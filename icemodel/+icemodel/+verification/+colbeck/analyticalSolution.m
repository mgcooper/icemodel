function solution = analyticalSolution(experiment_name, def)
   %ANALYTICALSOLUTION Compute the analytical Colbeck infiltration solution.
   %
   %  solution = icemodel.verification.colbeck.analyticalSolution(experiment)
   %  solution = icemodel.verification.colbeck.analyticalSolution(experiment, def)
   %
   %  Returns the closed-form column-integrated liquid water storage and bottom
   %  outflow for one of the three Colbeck 1976 / Clark 2017 infiltration
   %  experiments. The same kernel parameters as icemodel.column.liquid_flux are
   %  used (m_exp from parameterLookup) so analytical and numerical IceModel
   %  candidates share one constitutive law.
   %
   %  Methods
   %    Ripe snow (exp1):
   %       Isothermal kinematic-wave solution. A Rankine-Hugoniot shock advances
   %       at speed c_shock = (q_top - q(f_liq_0)) / (f_liq_inf - f_liq_0) until
   %       it reaches the column bottom; then storage plateaus at f_liq_inf *
   %       total_depth until rain shuts off. After shutoff a centered
   %       rarefaction at the top carries monotone-decreasing saturations to the
   %       bottom.
   %
   %    Cold snow (exp2 and exp3):
   %       Clark 2017 method-of-characteristics solution (Eqs. 10-13).
   %       The wetting front advances at
   %          dz_w/dt = q_top / (f_frz + f_liq_w - f_liq_0)
   %       where f_frz = (ro_ice * cp_ice / (ro_liq * Lf)) * f_ice * (Tf - T) is
   %       the thermal water requirement (Clark Eq. 10) and f_liq_w is the
   %       steady-state liquid fraction in the wetted region. Storage grows
   %       linearly with the front. After the front reaches the column bottom
   %       (or rain shuts off), recession follows a centered rarefaction; for
   %       partial-duration cases (exp3) this includes the merged shock per
   %       Clark Eqs. 28-32.
   %
   %  Inputs
   %    experiment_name  "exp1" / "exp2" / "exp3" matching caseDefinition.
   %    def              Optional pre-loaded case definition.
   %
   %  Outputs
   %    solution   Struct with fields:
   %       time_seconds                   n x 1 sample times [s]
   %       time_datetime                  n x 1 datetime stamps
   %       snow_liquid_water_storage_m    column-integrated liquid [m]
   %       bottom_outflow_mps             outflow at z = total_depth [m s-1]
   %       parameters                     struct of derived parameters
   %
   %  References
   %    Colbeck, S. C. (1976). An analysis of water flow in dry snow.
   %    Clark, M. P. et al. (2017). An analytical test case for snow models.
   %    Clark, M. P. et al. (2021). The numerical implementation of land
   %      models: Problem formulation and laugh tests.
   %
   % See also: icemodel.column.liquid_flux,
   %  icemodel.column.saturated_hydraulic_conductivity,
   %  icemodel.parameterLookup, icemodel.physicalConstant

   arguments
      experiment_name (1, 1) string
      def = icemodel.verification.colbeck.caseDefinition()
   end

   names = string({def.experiments.name});
   idx = find(names == experiment_name, 1);
   if isempty(idx)
      error('icemodel:colbeck:unknownExperiment', ...
         'unknown Colbeck experiment: %s', experiment_name);
   end
   experiment = def.experiments(idx);

   n_steps = def.time.n_steps;
   dt = def.time.dt_s;
   time_seconds = (1:n_steps).' * dt;
   time_datetime = datetime(1990, 1, 1, 0, 0, 0, 'TimeZone', 'UTC') ...
      + seconds(time_seconds);

   [m_exp, f_res_pore] = icemodel.parameterLookup('m_exp', 'f_res_pore');
   [Tf, Lf, ro_ice, ro_liq, cp_ice] = icemodel.physicalConstant( ...
      'Tf', 'Lf', 'ro_ice', 'ro_liq', 'cp_ice');

   q_top = def.boundary.q_top_m_per_s;
   total_depth = def.grid.total_depth_m;
   rain_window = def.time.rain_window_s;
   f_ice_0 = experiment.f_ice;
   f_liq_0 = experiment.f_liq;
   T_0 = experiment.T_K;

   % Use the case-specific intrinsic permeability (matches the SUMMA Colbeck
   % parameter trial value), so analytical and numerical IceModel candidates and
   % the SUMMA reference share one constitutive constant. The "darcy" k_sat
   % method is used throughout this file: k_sat = permeability * ro_liq * g / mu
   permeability = NaN;
   if isfield(experiment, 'permeability_m2')
      permeability = experiment.permeability_m2;
   end

   if f_liq_0 > 0
      % Ripe snow (exp1): isothermal kinematic-wave shock + recession.
      [storage, q_bot, params] = ripe_solution( ...
         time_seconds, q_top, rain_window, total_depth, ...
         f_ice_0, f_liq_0, f_res_pore, m_exp, permeability);
   else
      % Cold cases: closed-form Clark 2017 wetting-front solution covers
      % full-duration regimes (front reaches bottom before rain ends).
      % Partial-duration cases have a merged-shock structure (Clark Eqs. 28-32)
      % that does not admit a clean closed form; we solve the kinematic-wave PDE
      % on a fine grid as the analytical reference. Auto-route based on
      % t_arrival vs rain_window.
      [storage, q_bot, params] = cold_solution( ...
         time_seconds, q_top, rain_window, total_depth, f_ice_0, f_liq_0, ...
         T_0, f_res_pore, m_exp, Tf, Lf, ro_ice, ro_liq, cp_ice, permeability);

      if isfield(params, 'kind') && params.kind == "cold" ...
            && params.t_arrival_s > rain_window

         % Partial-duration: use fine-grid kinematic-wave reference.
         [storage, q_bot] = cold_partial_duration_pde( ...
            time_seconds, q_top, rain_window, total_depth, f_ice_0, T_0, ...
            f_res_pore, m_exp, Tf, Lf, ro_ice, ro_liq, cp_ice, permeability);
         params.kind = "cold_partial_pde";
      end
   end

   solution.time_seconds = time_seconds;
   solution.time_datetime = time_datetime;
   solution.snow_liquid_water_storage_m = storage;
   solution.bottom_outflow_mps = q_bot;
   solution.parameters = params;
end

%% Ripe snow: Colbeck 1972 / Clark 2017 Eqs. 5-9 closed-form solution.
function [storage, q_bot, params] = ripe_solution( ...
      time_seconds, q_top, rain_window, total_depth, ...
      f_ice_0, f_liq_0, f_res_pore, m_exp, permeability)

   % Pore / residual / capacity for the ripe column (uniform initial state, no
   % refreezing).
   f_por = 1 - f_ice_0;
   f_res = f_res_pore * f_por;
   avail = max(eps, f_por - f_res);

   % Saturated hydraulic conductivity from the darcy form using the
   % case-specific permeability. Returns column-shaped vector; pick the
   % (uniform) scalar value here for the closed-form derivation.
   k_sat = icemodel.column.saturated_hydraulic_conductivity( ...
      f_ice_0, f_liq_0, method="darcy", permeability=permeability);
   k_sat = k_sat(1);

   % S_inf is the relative saturation at infinity (the asymptotic equilibrium
   % saturation that balances q_top through the constitutive law). f_liq_inf and
   % c_inf below carry the same "at infinity / asymptotic equilibrium" meaning.
   S_inf = (q_top / k_sat) ^ (1 / m_exp);
   f_liq_inf = f_res + S_inf * avail;

   % Rankine-Hugoniot shock speed; q at f_liq_0 is zero by design (residual).
   q_at_0 = q_at_f_liq(f_liq_0, k_sat, f_res, avail, m_exp);
   c_shock = (q_top - q_at_0) / max(f_liq_inf - f_liq_0, eps);

   % Time when the shock reaches z = total_depth, and when the leading edge of
   % the post-rain rarefaction reaches the bottom.
   t_arrival = total_depth / c_shock;
   c_inf = m_exp * k_sat * S_inf ^ (m_exp - 1) / avail;
   t_drain_start = rain_window + total_depth / c_inf;

   n = numel(time_seconds);
   storage = zeros(n, 1);
   q_bot = zeros(n, 1);
   for k = 1:n
      [storage(k), q_bot(k)] = sample_ripe(time_seconds(k), q_top, ...
         rain_window, t_arrival, t_drain_start, f_liq_0, f_liq_inf, ...
         total_depth, c_inf, k_sat, avail, m_exp);
   end

   params = struct( ...
      'kind',           "ripe", ...
      'k_sat_m_per_s',  k_sat, ...
      'f_por',          f_por, ...
      'f_res',          f_res, ...
      'f_liq_0',        f_liq_0, ...
      'f_liq_inf',      f_liq_inf, ...
      'S_inf',          S_inf, ...
      'c_shock_m_per_s',c_shock, ...
      'c_inf_m_per_s',  c_inf, ...
      't_arrival_s',    t_arrival, ...
      't_drain_start_s',t_drain_start);
end

function [storage, q_bot] = sample_ripe(t, q_top, rain_window, ...
      t_arrival, t_drain_start, f_liq_0, f_liq_inf, total_depth, ...
      c_inf, k_sat, avail, m_exp)

   if t <= 0
      storage = f_liq_0 * total_depth;
      q_bot = 0;
      return
   end

   if t <= rain_window
      if t < t_arrival
         z_front = total_depth * (t / t_arrival);
         storage = f_liq_0 * total_depth + (f_liq_inf - f_liq_0) * z_front;
         q_bot = 0;
      else
         storage = f_liq_inf * total_depth;
         q_bot = q_top;
      end
      return
   end

   if t < t_drain_start
      tau = t - rain_window;
      storage = f_liq_inf * total_depth - q_top * tau;
      q_bot = q_top;
      return
   end

   % Recession via the characteristic that arrives at the bottom at t.
   c_at_bot = total_depth / (t - rain_window);
   if c_at_bot >= c_inf
      S_star = (q_top / k_sat) ^ (1 / m_exp);
   else
      S_star = max(0, c_at_bot * avail / (m_exp * k_sat)) ^ (1 / (m_exp - 1));
   end
   q_bot = k_sat * S_star ^ m_exp;

   % Storage from analytical mass-balance integral; closed form for m=2,
   % numerical quadrature otherwise.
   if abs(m_exp - 2) < eps
      pref = (avail ^ 2) / (4 * k_sat);
      D2 = total_depth ^ 2;
      tau1 = t_drain_start - rain_window;
      tau = t - rain_window;
      drained = pref * D2 * (1 / tau1 - 1 / tau);
   else
      n = 256;
      s = linspace(t_drain_start, t, n);
      qs = zeros(size(s));
      for i = 1:n
         c_at_bot_s = total_depth / (s(i) - rain_window);
         S_s = max(0, c_at_bot_s * avail / (m_exp * k_sat)) ^ (1 / (m_exp - 1));
         qs(i) = k_sat * S_s ^ m_exp;
      end
      drained = trapz(s, qs);
   end
   storage_at_drain_start = f_liq_inf * total_depth ...
      - q_top * (t_drain_start - rain_window);
   storage = max(0, storage_at_drain_start - drained);
end

%% Cold snow: Clark 2017 wetting-front advance with phase-change uptake.
function [storage, q_bot, params] = cold_solution( ...
      time_seconds, q_top, rain_window, total_depth, ...
      f_ice_0, f_liq_0, T_0, f_res_pore, m_exp, ...
      Tf, Lf, ro_ice, ro_liq, cp_ice, permeability)

   % Thermal water requirement (Clark 2017 Eq. 10):
   %   f_frz = (ro_ice * cp_ice / (ro_liq * Lf)) * f_ice * (Tf - T)
   % This is the volumetric liquid that must freeze to bring the dry layer to
   % Tf. After freezing, the layer is at Tf and behaves as ripe snow.
   f_frz = (ro_ice * cp_ice / (ro_liq * Lf)) * f_ice_0 * (Tf - T_0);

   % After thermal saturation the wetted region carries the steady- state
   % saturation that balances q_top through the post-refreeze ice fraction.
   % Refreezing increases f_ice by
   %   f_frz / (ro_ice / ro_liq) = f_frz * ro_liq / ro_ice.
   f_ice_w = f_ice_0 + f_frz * (ro_liq / ro_ice);
   f_por_w = 1 - f_ice_w;
   f_res_w = f_res_pore * f_por_w;
   avail_w = max(eps, f_por_w - f_res_w);

   % Saturated hydraulic conductivity in the wetted region (darcy form with
   % case-specific permeability).
   k_sat_w_vec = icemodel.column.saturated_hydraulic_conductivity( ...
      f_ice_w, 0, method="darcy", permeability=permeability);
   k_sat_w = k_sat_w_vec(1);
   S_w = min(1, (q_top / k_sat_w) ^ (1 / m_exp));
   f_liq_w = f_res_w + S_w * avail_w;

   % Wetting-front advance speed (Clark 2017 Eq. 13):
   %   dz_w/dt = q_top / (f_frz + f_liq_w - f_liq_0)
   uptake = max(eps, f_frz + f_liq_w - f_liq_0);
   c_w = q_top / uptake;
   t_arrival = total_depth / c_w;

   % Total inflow at the top, by case (used for phase 1b mass balance).
   inflow_total_at_rain_end = q_top * rain_window;

   % Recession parameters once the wetted region is mature.
   c_inf = m_exp * k_sat_w * S_w ^ (m_exp - 1) / avail_w;
   t_drain_start = min(rain_window, t_arrival) + total_depth / c_inf;

   n = numel(time_seconds);
   storage = zeros(n, 1);
   q_bot = zeros(n, 1);
   for k = 1:n
      [storage(k), q_bot(k)] = sample_cold(time_seconds(k), q_top, ...
         rain_window, t_arrival, t_drain_start, f_liq_0, f_liq_w, ...
         total_depth, c_w, c_inf, k_sat_w, f_res_w, avail_w, m_exp, ...
         f_frz, inflow_total_at_rain_end);
   end

   params = struct( ...
      'kind',            "cold", ...
      'k_sat_m_per_s',   k_sat_w, ...
      'f_frz',           f_frz, ...
      'f_ice_w',         f_ice_w, ...
      'f_por_w',         f_por_w, ...
      'f_res_w',         f_res_w, ...
      'f_liq_w',         f_liq_w, ...
      'S_w',             S_w, ...
      'c_w_m_per_s',     c_w, ...
      'c_inf_m_per_s',   c_inf, ...
      't_arrival_s',     t_arrival, ...
      't_drain_start_s', t_drain_start);
end

function [storage, q_bot] = sample_cold(t, q_top, rain_window, ...
      t_arrival, t_drain_start, f_liq_0, f_liq_w, total_depth, ...
      c_w, c_inf, k_sat_w, f_res, avail, m_exp, f_frz, ...
      inflow_total_at_rain_end) %#ok<INUSL>

   if t <= 0
      storage = f_liq_0 * total_depth;
      q_bot = 0;
      return
   end

   % Effective time at which the column becomes fully wet at f_liq_w. For
   % full-duration cases (front reaches bottom before rain ends), this is
   % t_arrival; for partial-duration cases (rain ends before arrival), gravity
   % drainage continues to feed the front advance until it reaches the bottom,
   % so the column reaches a fully wet state at t_arrival.
   t_full = t_arrival;

   % Phase 1: wetting front still advancing. No bottom outflow until the front
   % reaches z = total_depth.
   %
   % Phase 1a (during rain, t <= rain_window): the wetted region behind the
   % shock is at f_liq_w by Rankine-Hugoniot mass balance. Storage therefore
   % equals f_liq_w * z_w(t).
   %
   % Phase 1b (rain stopped but front still advancing): no inflow at the top, so
   % storage drops via the thermal-saturation cost as the front advances. By
   % mass balance,
   %   storage(t) = inflow_total_at_rain_end - f_frz * z_w(t)
   % The wet region contracts at the top via a rarefaction even as the front
   % advances at the bottom (Clark 2017 Eqs. 28-32 merged shock); this
   % mass-balance form captures the integrated effect without requiring the full
   % closed-form merged-shock advance.
   if t < t_full
      z_w = min(total_depth, c_w * t);
      if t <= rain_window
         storage = f_liq_w * z_w + f_liq_0 * (total_depth - z_w);
      else
         storage = max(0, inflow_total_at_rain_end - f_frz * z_w);
      end
      q_bot = 0;
      return
   end

   % Phase 2: column fully wet at f_liq_w. Recession follows once the leading
   % rarefaction from the top reaches the bottom.
   tau = t - t_full;
   if rain_window > t_full
      % Full-duration: rain still active beyond t_full.
      tau_after_rain = max(0, t - rain_window);
   else
      % Partial-duration: rain stopped at rain_window before t_full.
      % Effective rarefaction starts at t_full.
      tau_after_rain = tau;
   end

   t_off_recession = max(rain_window, t_full);
   % For partial-duration cases (rain ended before column was fully wetted), the
   % rarefaction at the top has already reached the bottom by the time the
   % wetting front arrives, so there is no f_liq_w plateau phase. Recession
   % starts immediately at t_full.
   if t_full <= rain_window
      t_drain_start_eff = t_off_recession + total_depth / c_inf;
   else
      t_drain_start_eff = t_off_recession;
   end

   if t < rain_window
      % Should not occur given Phase 1 covers t < t_full and
      % t_full = t_arrival >= rain_window or t_full < rain_window.
      storage = f_liq_w * total_depth;
      q_bot = q_top;
      return
   end

   % Storage at the moment the column is fully wetted (t = t_full). For
   % full-duration cases (t_full < rain_window) the column is at f_liq_w
   % throughout. For partial-duration the wet region above has rarefied during
   % phase 1b; mass balance gives:
   %   storage_at_full = inflow_total_at_rain_end - f_frz * total_depth
   % which is also equal to f_liq_w * total_depth when phase 1b is empty
   % (full-duration case), since then inflow = q_top * t_arrival and q_top = c_w
   % * (f_liq_w + f_frz), giving f_liq_w * D as expected.
   if t_full <= rain_window
      storage_at_full = f_liq_w * total_depth;
   else
      storage_at_full = max(0, inflow_total_at_rain_end - f_frz * total_depth);
   end

   if t < t_drain_start_eff
      % Storage drains at q_top (full-duration: rain still on or just off;
      % outflow tracks q_top). For partial-duration, gravity drainage at q_w
      % continues and outflow is q_top until the rarefaction reaches the bottom.
      storage = max(0, storage_at_full - q_top * tau_after_rain);
      q_bot = q_top;
      return
   end

   % Recession via the characteristic that arrives at the bottom at t.
   c_at_bot = total_depth / (t - t_off_recession);
   if c_at_bot >= c_inf
      S_star = (q_top / k_sat_w) ^ (1 / m_exp);
   else
      S_star = max(0, c_at_bot * avail / (m_exp * k_sat_w)) ^ (1 / (m_exp - 1));
   end
   q_bot = k_sat_w * S_star ^ m_exp;

   % Storage via mass-balance integral (numerical quadrature).
   drained = drained_quad(t_off_recession, t, t_drain_start_eff, ...
      total_depth, k_sat_w, avail, m_exp);
   storage_at_drain_start = storage_at_full ...
      - q_top * (t_drain_start_eff - max(rain_window, t_full));
   storage = max(0, storage_at_drain_start - drained);
end

function drained = drained_quad(t_off, t_now, t_drain_start, ...
      total_depth, k_sat, avail, m_exp)
   n = 256;
   s = linspace(t_drain_start, t_now, n);
   qs = zeros(size(s));
   for i = 1:n
      c_at_bot = total_depth / (s(i) - t_off);
      S_star = max(0, c_at_bot * avail / (m_exp * k_sat)) ^ (1 / (m_exp - 1));
      qs(i) = k_sat * S_star ^ m_exp;
   end
   drained = trapz(s, qs);
end

function q = q_at_f_liq(f_liq, k_sat, f_res, avail, m_exp)
   if avail <= 0 || f_liq <= f_res
      q = 0;
      return
   end
   S = (f_liq - f_res) / avail;
   q = k_sat * S ^ m_exp;
end

%% Cold partial-duration: kinematic-wave PDE on a fine grid.
%
% Solves the same constitutive law as the production kernel (q = k_sat * S^m_exp
% with thermal-saturation uptake at the wetting front) using an explicit upwind
% scheme on a 10x-refined grid. This is the analytical reference for cases where
% the closed-form merged-shock solution (Clark 2017 Eqs. 28-32) is awkward to
% derive by hand: it solves the same PDE the closed form would, just numerically
% and at higher spatial / temporal resolution than the production candidate.
%
% Independence from the production kernel: this function uses only explicit FV
% upwind on f_liq with an in-place refreezing step against local cold content,
% so changes to icemodel.column.infiltration cannot mask issues caught here.

function [storage, q_bot] = cold_partial_duration_pde( ...
      time_seconds, q_top, rain_window, total_depth, f_ice_0, T_0, ...
      f_res_pore, m_exp, Tf, Lf, ro_ice, ro_liq, cp_ice, permeability)

   % Refined grid and timestep. 1000 layers x 1-second substep gives 10x grid
   % refinement and 60x time refinement vs the production candidate (100 layers,
   % ~7 s substep on the Colbeck case).
   N_ref  = 1000;
   dz_ref = total_depth / N_ref;
   dt_ref = 1.0;

   % Initial column state (uniform; cold-snow case starts dry).
   T     = T_0     * ones(N_ref, 1);
   f_ice = f_ice_0 * ones(N_ref, 1);
   f_liq = zeros(N_ref, 1);

   % Constant k_sat from the darcy form (precomputed once because method =
   % "darcy" returns a uniform value driven by permeability).
   k_sat_vec = icemodel.column.saturated_hydraulic_conductivity( ...
      f_ice, f_liq, method="darcy", permeability=permeability);
   k_sat_layers = k_sat_vec(1) * ones(N_ref, 1);

   % Ice-to-liquid water-equivalent ratio; needed for the f_liq upper-bound
   % clip after each substep.
   ro_iwe = ro_ice / ro_liq;

   % Per-output-step diagnostics. dt_main is the requested output cadence (60 s
   % for Colbeck); each main step takes n_sub fine substeps of dt_ref.
   n_steps = numel(time_seconds);
   dt_main = time_seconds(2) - time_seconds(1);
   storage = zeros(n_steps, 1);
   q_bot   = zeros(n_steps, 1);

   for step = 1:n_steps

      % Top boundary follows the Colbeck rain window.
      t_end = time_seconds(step);
      q_top_now = q_top * (t_end <= rain_window + eps);

      n_sub = round(dt_main / dt_ref);
      bot_flux_accum = 0;

      for sub = 1:n_sub

         % --- Pore / residual / capacity / relative saturation -------
         f_por    = 1 - f_ice;
         f_res    = f_res_pore .* f_por;
         availCap = max(eps, f_por - f_res);
         relSat   = max(0, (f_liq - f_res) ./ availCap);
         relSat   = min(1, relSat);

         % --- Per-layer outgoing flux (upwind) -----------------------
         %   q_layer(i) = flux out the bottom of layer i, computed at
         %               layer i's own state (S, k_sat).
         q_layer = k_sat_layers .* relSat .^ m_exp;

         % --- Boundary conditions / interface flux vector -----------
         %   bc_N = q_top_now    (surface inflow)
         %   bc_S = q_layer(end) (free-drainage at the bottom)
         q_int = [q_top_now; q_layer];
         q_bot_sub = q_layer(end);

         % --- Net liquid water flux per layer ------------------------
         %   d(f_liq) = (dt_ref / dz_ref) * d/dz q(z) [1]
         dliq = (dt_ref / dz_ref) * (q_int(1:end-1) - q_int(2:end));
         f_liq = f_liq + dliq;

         % --- Clip f_liq to physical bounds --------------------------
         %   Lower: f_liq >= 0 (no negative liquid).
         %   Upper: f_liq <= ro_iwe * (1 - f_ice) (cannot exceed pore
         %          volume in liquid-water-equivalent units).
         f_liq = max(0, f_liq);
         f_air = ro_iwe * (1 - f_ice) - f_liq;
         neg = f_air < 0;
         if any(neg)
            f_liq(neg) = ro_iwe * (1 - f_ice(neg));
         end

         bot_flux_accum = bot_flux_accum + q_bot_sub * dt_ref;

         % --- Cold-content refreezing per layer ----------------------
         % Per layer with T < Tf and f_liq > 0, freeze up to the available cold
         % content, releasing latent heat to warm the layer toward Tf. Same
         % physics as icemodel.column.infiltration's refreezing block.
         for k = 1:N_ref
            if T(k) >= Tf - eps || f_liq(k) <= 0 || f_ice(k) <= 0
               continue
            end
            cc = ro_ice * f_ice(k) * cp_ice * (Tf - T(k));
            df_frz = min(f_liq(k), cc / (ro_liq * Lf));
            if df_frz <= 0
               continue
            end
            f_liq(k) = f_liq(k) - df_frz;
            f_ice(k) = f_ice(k) + df_frz / ro_iwe;
            dT = (ro_liq * df_frz * Lf) / (ro_ice * f_ice(k) * cp_ice);
            T(k) = min(Tf, T(k) + dT);
         end
      end

      % Step diagnostics (column-integrated liquid storage and step- averaged
      % bottom outflow).
      storage(step) = sum(f_liq) * dz_ref;
      q_bot(step)   = bot_flux_accum / dt_main;
   end
end
