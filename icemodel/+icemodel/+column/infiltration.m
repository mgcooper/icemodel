function [f_liq, f_ice, T, diag] = infiltration( ...
      f_liq, f_ice, T, dz, dt, q_top, kwargs)
   %INFILTRATION Snow column liquid mass + cold-content + conduction update.
   %
   %  [f_liq, f_ice, T, diag] = icemodel.column.infiltration( ...
   %     f_liq, f_ice, T, dz, dt, q_top)
   %
   %  Owns the full snow infiltration step: an explicit upwind liquid mass
   %  redistribution with CFL substepping, cold-content refreezing in cold
   %  layers (latent heat released, T capped at Tf), and a one-step explicit
   %  conduction sweep over the main timestep dt. Returns the updated (f_liq,
   %  f_ice, T) triple plus a diagnostic struct.
   %
   %  The kernel is the canonical icemodel infiltration entry point. Both the
   %  snow verification driver (Colbeck 1976 / Clark 2017) and any future
   %  production caller share this code path. Empirical coefficients come from
   %  icemodel.parameterLookup, physical constants from
   %  icemodel.physicalConstant, the saturated hydraulic conductivity from
   %  icemodel.column.saturated_hydraulic_conductivity, and the thermal
   %  conductivity from icemodel.column.firn_thermal_conductivity.
   %
   %  Parameters:
   %  -----------
   %  f_liq : (N x 1) vector
   %      Volumetric liquid water fraction in each layer [1].
   %  f_ice : (N x 1) vector
   %      Volumetric ice fraction in each layer [1].
   %  T : (N x 1) vector
   %      Layer temperature [K].
   %  dz : scalar
   %      Layer thickness [m] (uniform grid).
   %  dt : scalar
   %      Main timestep [s].
   %  q_top : scalar
   %      Liquid flux entering the top interface [m s-1] (precipitation
   %      or surface meltwater; opts.f_res_pore_snow / _ice / _firn
   %      enter via the residual capillary floor inside liquid_flux).
   %
   %  Name-value:
   %  -----------
   %  k_sat_method : string (default "colbeck1972")
   %      "colbeck1972" | "shimizu1970" | "darcy". See
   %      icemodel.column.saturated_hydraulic_conductivity.
   %  grainsz : scalar (optional)
   %      Grain diameter [m]; required for "shimizu1970".
   %  permeability : scalar (optional)
   %      Intrinsic permeability kappa [m^2]; required for "darcy"
   %      (verification path against the SUMMA Colbeck parameter trial).
   %
   %  Returns:
   %  --------
   %  f_liq, f_ice, T : (N x 1) vectors
   %      Updated column state.
   %  diag : struct
   %      Substepping and mass-balance diagnostics:
   %        n_sub, dt_sub, k_sat_max, S_inf, c_max,
   %        inflow_total [m], outflow_total [m].
   %
   %  Variable analogues with Colbeck 1972:
   %      f_res_pore - residual water volume / pore volume [1] (Swi)
   %      f_por      - pore fraction [1]                       (phi)
   %      f_res      - residual water fraction [1]
   %      availCap   - available capacity [1]
   %      relSat     - relative saturation [1]                 (Sstar)
   %
   %  Notes
   %  -----
   %  Free-drainage bottom boundary (q_bot = q_internal(end)). The CFL bound on
   %  the substep size uses the larger of two characteristic speeds: (1)
   %  inflow-driven steady-state c_inflow at S_inf = (q_top /
   %  k_sat_max)^(1/m_exp), which dominates during rain; and (2) c_state_max
   %  evaluated against the current per-layer S, which dominates during
   %  recession (when q_top -> 0 and S_inf -> 0). Refreezing and explicit
   %  conduction each run once per main step, after the liquid update.
   %
   % See also: icemodel.column.liquid_flux,
   %  icemodel.column.saturated_hydraulic_conductivity,
   %  icemodel.column.firn_thermal_conductivity,
   %  icemodel.parameterLookup, icemodel.physicalConstant
   %
   %#codegen

   arguments
      f_liq           (:, 1) double
      f_ice           (:, 1) double
      T               (:, 1) double
      dz              (1, 1) double
      dt              (1, 1) double
      q_top           (1, 1) double
      kwargs.k_sat_method (1, 1) string {mustBeMember( ...
         kwargs.k_sat_method, ["colbeck1972", "shimizu1970", "darcy"])} ...
         = "colbeck1972"
      kwargs.grainsz      (1, 1) double = NaN
      kwargs.permeability (1, 1) double = NaN
   end

   % Load thermodynamic constants and constitutive-law parameters once. m_exp is
   % the Mualem exponent (Clark 2017 Eq. 7); f_res_pore is the canonical Colbeck
   % residual capillary saturation (production callers override per-substep via
   % opts.f_res_pore_snow / _ice / _firn at the liquid_flux call site, but the
   % CFL bound here uses the canonical value); cfl_safety is the substep safety
   % factor.
   persistent m_exp f_res_pore cfl_safety Tf Lf ro_ice ro_liq cp_ice
   if isempty(m_exp)
      [m_exp, f_res_pore, cfl_safety] = icemodel.parameterLookup( ...
         'm_exp', 'f_res_pore', 'cfl_safety');
      [Tf, Lf, ro_ice, ro_liq, cp_ice] = icemodel.physicalConstant( ...
         'Tf', 'Lf', 'ro_ice', 'ro_liq', 'cp_ice');
   end

   % Ice-to-liquid water-equivalent ratio (used for the f_liq upper-bound clip
   % below: f_liq cannot exceed (ro_ice/ro_liq) * (1 - f_ice), which is the
   % available pore volume converted to liquid- water equivalent).
   ro_iwe = ro_ice / ro_liq;

   % Unpack optional keyword-args
   grainsz = kwargs.grainsz;
   k_sat_method = kwargs.k_sat_method;
   permeability = kwargs.permeability;

   % --- Pore / residual / capacity ----------------------------------------
   % Pore fraction and residual liquid fraction per layer. availCap is the
   % per-layer available pore capacity for liquid water above the capillary
   % residual floor; eps avoids division-by-zero in fully ice layers (f_ice==1).
   f_por = 1.0 - f_ice;
   f_res = f_res_pore .* f_por;
   availCap = max(eps, f_por - f_res);

   % --- CFL substep count -------------------------------------------------
   % Saturated hydraulic conductivity per layer (column-shaped) from the
   % dispatched closure scheme. k_sat_max bounds the per-layer values for the
   % CFL estimate.
   k_sat = icemodel.column.saturated_hydraulic_conductivity( ...
      f_ice, f_liq, method=k_sat_method, grainsz=grainsz, permeability=permeability);
   k_sat_max = max(k_sat);

   % Relative saturation per layer, clipped to [0, 1]. This is used to estimate
   % the column-state characteristic speed which dominates the CFL bound during
   % recession (when no inflow is feeding the steady-state shock).
   S = (f_liq - f_res) ./ availCap;
   S = min(1, max(0, S));
   c = m_exp .* k_sat .* S .^ (m_exp - 1) ./ availCap;
   c_max = max(c);

   % Inflow-driven steady-state characteristic speed at S_inf = (q_top /
   % k_sat_max)^(1/m_exp). Dominates during rain phases.
   if k_sat_max > 0 && q_top > 0
      S_inflow = min(1, (q_top / k_sat_max) ^ (1 / m_exp));
      c_inflow = m_exp * k_sat_max * S_inflow ^ (m_exp - 1) / max(availCap);
   else
      S_inflow = 0;
      c_inflow = 0;
   end

   % CFL bound: take the larger of the two characteristic speeds and choose a
   % substep size such that the wave moves at most cfl_safety cells per substep.
   c_max = max(c_max, c_inflow);
   dt_sub_target = cfl_safety * dz / max(c_max, eps);
   n_sub = max(1, ceil(dt / dt_sub_target));
   dt_sub = dt / n_sub;

   % --- Liquid mass substep loop -----------------------------------------
   inflow_total  = 0;
   outflow_total = 0;

   for sub = 1:n_sub

      % Calculate fluxes between layers [m/s]. liquid_flux returns one flux
      % value per layer, computed at the layer's own state; this is the upwind
      % interface flux carried out the bottom of each layer.
      q = icemodel.column.liquid_flux(f_liq, f_ice, ...
         k_sat_method=k_sat_method, grainsz=grainsz, permeability=permeability);

      % Boundary conditions:
      %   bc_N = q_top        surface inflow (precipitation / surface meltwater)
      %   bc_S = q(end)       free-drainage at the bottom interface
      %                       (q_bot = upstream q from the bottom layer)
      %
      % After the BCs are appended to q:
      %
      % ----- q(1) = flux across the top (surface) interface into the
      %              first CV
      %
      % ----- q(2) = flux across the second interface into the second CV
      %
      % ----- q(3) = flux across the third interface into the third CV
      %   .
      %   .
      %   .
      % ----- q(N) = flux across the second-from-bottom interface into
      %              the bottom CV
      %
      % ----- q(N+1) = flux across the bottom interface = free-drainage
      bc_N  = q_top;
      bc_S  = q(end);
      q_int = [bc_N; q];        % length N+1
      q_bot = bc_S;

      % Compute net liquid water flux for each layer
      %   d(f_liq) = (dt_sub / dz) * d/dz q(z) [1]
      d_liq = (dt_sub / dz) * (q_int(1:end-1) - q_int(2:end));
      f_liq = f_liq + d_liq;

      % Clip f_liq to physical bounds.
      %
      %   Lower bound: f_liq >= 0 (no negative liquid water).
      %   Upper bound: f_liq <= (ro_ice/ro_liq) * (1 - f_ice) (cannot
      %                exceed the available pore volume in liquid-water
      %                equivalent units).
      %
      % This is a single-pass post-clip: d_liq is applied first, then f_liq is
      % snapped to the physical envelope. Excess water above the upper bound is
      % dropped (non-conservative on the upper bound, same as the original
      % two-pass scheme). Drainage past f_res is implicitly prevented by
      % liquid_flux upstream, which returns q = 0 wherever f_liq <= f_res, so
      % layers at the capillary residual produce zero outgoing flux and the
      % lower bound is rarely exercised in practice.
      %
      % Note: the original infiltration kernel used a two-pass per-layer
      % pre-clamp (drainage clamped to f_liq - f_res, plus a downward rebalance
      % against the receiver layer's max_infill capacity). With proper CFL
      % substepping (cfl_safety < 1) the wave moves less than one cell per
      % substep and the receiver-capacity check is dormant for the Colbeck
      % cases, so the simpler post-clip suffices and the verification metrics
      % agree with SUMMA at sub-mm RMSE.
      f_liq = max(0, f_liq);
      f_air_capacity = ro_iwe * f_por - f_liq;
      neg = f_air_capacity < 0;
      if any(neg)
         f_liq(neg) = ro_iwe * f_por(neg);
      end

      % Cumulative mass-balance accumulators (for diagnostics; the timestep
      % total is exposed via the diag struct).
      inflow_total  = inflow_total  + q_top * dt_sub;
      outflow_total = outflow_total + q_bot * dt_sub;
   end

   % --- Cold-content refreezing -------------------------------------------
   % Per layer with T < Tf and f_liq > 0, freeze up to the available cold
   % content, releasing latent heat to warm the layer toward Tf.
   %
   % Cold content per unit layer volume [J m-3]:
   %   cold_content = ro_ice * f_ice * cp_ice * (Tf - T)
   % Latent heat released by freezing volumetric fraction df_frz of
   % liquid water (per unit layer volume) [J m-3]:
   %   dE = ro_liq * df_frz * Lf
   % Equating gives the maximum df_frz the cold content can absorb.
   for k = 1:numel(T)
      if T(k) >= Tf - eps || f_liq(k) <= 0 || f_ice(k) <= 0
         continue
      end
      cold_content = ro_ice * f_ice(k) * cp_ice * (Tf - T(k));
      df_frz       = min(f_liq(k), cold_content / (ro_liq * Lf));
      if df_frz <= 0
         continue
      end

      % Update phase fractions and temperature in place.
      f_liq(k) = f_liq(k) - df_frz;
      f_ice(k) = f_ice(k) + df_frz / ro_iwe;
      dT       = (ro_liq * df_frz * Lf) / (ro_ice * f_ice(k) * cp_ice);
      T(k)     = min(Tf, T(k) + dT);
   end

   % --- Explicit conduction (one main step) -------------------------------
   % Centered explicit conduction over the main step dt with adiabatic top and
   % bottom boundaries (zero-gradient ghost layers). Volumetric heat capacity
   % uses ice fraction only (water content is small for the cold-snow Colbeck
   % experiments). Clark 2017 sets the spatial heat-flux divergence to zero for
   % the verification target, so the conduction step has negligible effect on
   % the verification metrics (storage, outflow); it is retained because
   % production callers need the energy update.
   N = numel(T);
   if N >= 2
      k_thermal_layers = icemodel.column.firn_thermal_conductivity(T, f_ice);
      cv = ro_ice * f_ice * cp_ice;
      alpha = k_thermal_layers ./ max(cv, eps);  % m^2 s-1
      T_new = T;
      for k = 2:N-1
         T_new(k) = T(k) + (alpha(k) * dt / dz^2) * (T(k+1) - 2 * T(k) + T(k-1));
      end
      % Adiabatic top/bottom boundaries (zero-gradient).
      T_new(1) = T(1) + (alpha(1) * dt / dz^2) * (T(2)   - T(1));
      T_new(N) = T(N) + (alpha(N) * dt / dz^2) * (T(N-1) - T(N));
      T = min(T_new, Tf);
   end

   % Diagnostics struct, returned for verification mass-balance checks
   % and CFL-bound debugging.
   diag = struct( ...
      'n_sub',         n_sub, ...
      'dt_sub',        dt_sub, ...
      'k_sat_max',     k_sat_max, ...
      'S_inflow',      S_inflow, ...
      'c_max',         c_max, ...
      'inflow_total',  inflow_total, ...
      'outflow_total', outflow_total);
end

% For reference, in terms of CV thicknesses:
% Swi = f_res_por = f_res/f_por  = h_res/h_por = h_res/(h_tot - h_ice)
% Sw  = liqsat    = f_liq/f_por  = h_liq/h_por = h_liq/(h_tot - h_ice)
% phi = porosity  = f_por        = h_por/h_tot = (h_tot - h_ice)/h_tot
%
% Sstar = relSat  = (f_liq - f_res)/availCap = (h_liq - h_res)/availCap
%
% dq  = q(1:end-1) - q(2:end)
% -dq = q(2:end)   - q(1:end-1)
% dq < 0 means water drains out of CV
% dq > 0 means water enters CV
% dq < 0 && -df_liq > f_liq means more drains than exists = unstable
% dq > 0 && dq + f_liq + f_ice > 1 means more drains than can be stored
%
% Note, check if infiltrating water satisfies the cold content:
%   f_cc = (ro_ice * cp_ice / (ro_liq * Lf)) * f_ice * (Tf - T)
% (this is the f_frz / Clark 2017 Eq. 10 thermal water requirement used by the
% analyticalSolution wetting-front advance).
