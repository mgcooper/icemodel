function [T, f_ice, f_liq, k_eff, ok, iter, a1, err] = ...
      solve_column_enthalpy(T_sfc, T, f_ice, f_liq, Fc, Fp, Sc, Sp, dz, ...
      delz, fn, dt, solver, tol, maxiter, ~, ~, ~, debug, varargin)
   %SOLVE_COLUMN_ENTHALPY Solve the column enthalpy balance.
   %
   % The signature keeps the alpha, use_aitken, and jumpmax inputs to match the
   % thermal-solver option list. Node-wise Aitken acceleration is off here.
   %
   % A trailing true selects the coupled vapor mode. The bulk conductivity
   % then leaves the vapor term out. The vapor energy travels on the two
   % face terms from icemodel.column.vapor_face_conductance, built from the
   % same face quantities as the vapor mass flux: a positive matrix part
   % and a deferred-correction source flux. The second trailing argument is
   % f_res_por, which the donor-phase predicate needs. Without them the
   % vapor term stays inside k_eff on the node-tangent slope, which is the
   % default.
   %
   %#codegen

   % The coupled vapor mode is opt-in and nothing sets it by default.
   use_coupled_vapor = nargin > 19 && varargin{1};

   % The residual pore fraction rides only with the coupled mode; the
   % default path never reads it.
   f_res_por = 0.0;
   if nargin > 20
      f_res_por = varargin{2};
   end

   % Update the water fraction
   f_wat = icemodel.column.water_fraction(f_ice, f_liq);

   % Update the melt-zone volumetric liquid fraction bounds at T=TL and T=TH.
   [f_liq_min, f_liq_max] = icemodel.column.meltzone_bounds(f_wat);

   % Compute vapor density [kg m-3]
   ro_vap = icemodel.vapor.saturation_vapor_density(T, f_liq);

   % Compute enthalpy [J m-3]
   H_old = icemodel.column.bulk_enthalpy(T, f_ice, f_liq, f_wat, ro_vap);

   % Store past values
   T_ice_old = T;
   f_liq_old = f_liq;

   % Initialize current values
   T_iter = T_ice_old + 2 * tol;

   % Initial past Picard iterates for Aitken-acceleration (disabled)
   % T_1 = nan(size(T));
   % T_2 = nan(size(T));

   % Iterate to solve the nonlinear heat equation
   ok = false;
   for iter = 0:maxiter-1

      % Update vapor density and derivative [kg m-3, kg m-3 K-1]
      [ro_vap, dro_vapdT] = icemodel.vapor.saturation_vapor_density( ...
         T, f_liq);

      % Update vapor thermal conductivity [W m-1 K-1]
      [k_vap, De] = icemodel.vapor.vapor_thermal_conductivity( ...
         T, f_liq, dro_vapdT);

      % Coupled mode moves the vapor term from the node conductivity to the
      % faces. The energy the solve carries is then the mass flux times the
      % donor cell's latent heat: the positive matrix part joins the
      % interface conductivity and the deferred flux joins the source, so
      % the converged face flux is exact.
      if use_coupled_vapor
         [k_vap_faces, q_vap_deferred] = ...
            icemodel.column.vapor_face_conductance( ...
            T, f_ice, f_liq, ro_vap, dro_vapdT, De, delz, fn, f_res_por);
         k_vap = zeros(size(k_vap));
      end

      % Update bulk thermal conductivity
      k_eff = icemodel.column.bulk_thermal_conductivity( ...
         T, f_ice, f_liq, k_vap);

      % Update bulk enthalpy and derivative wrt temperature
      [H, dHdT, dFdT] = icemodel.column.bulk_enthalpy( ...
         T, f_ice, f_liq, f_wat, ro_vap, dro_vapdT);

      % Update the general equation coefficients
      if use_coupled_vapor
         [aN, aP, aS, b, iM, a1, a2] = ...
            icemodel.column.assemble_enthalpy_system( ...
            T, f_ice, f_liq, dHdT, dFdT, dro_vapdT, H - H_old, Sc, Sp, ...
            k_eff, delz, fn, dz, dt, T_sfc, Fc, Fp, solver, ...
            k_vap_faces, q_vap_deferred);
      else
         [aN, aP, aS, b, iM, a1, a2] = ...
            icemodel.column.assemble_enthalpy_system( ...
            T, f_ice, f_liq, dHdT, dFdT, dro_vapdT, H - H_old, Sc, Sp, ...
            k_eff, delz, fn, dz, dt, T_sfc, Fc, Fp, solver);
      end

      % % Check diagonal dominance and condition number
      % icemodel.checkdiags(aP, aN, aS)

      % Exit here so the state variables are updated on the final iteration
      if all(abs(T - T_iter) < tol)
         ok = true;
         break
      end

      % Capture past values
      T_iter = T;

      % Solve the equation (predictor step)
      T = icemodel.numerics.trisolve(-aN, aP, -aS, b);

      % Update the temperature-enthalpy relationship (corrector step)
      [T, f_ice, f_liq, ok_mz] = icemodel.column.meltzone_transform( ...
         T, T_iter, f_liq, f_wat, dFdT, f_liq_min, f_liq_max, iM, debug);

      % If failure, return to the main program and shorten the timestep
      if ~ok_mz
         if debug
            dumpIceEnbalFailure("mztransform_rejected_state", T, ...
               T_ice_old, T_iter, f_ice, f_liq, f_wat, k_eff, Sc, dt, ...
               T_sfc, iM, iter, maxiter, aN, aP, aS, b);
         end
         return
      end

      % Control volume check - max water cannot exceed ro_ice / ro_liq.
      assertF(@() icemodel.column.assert_max_water(f_ice, f_liq));

      % Relaxation and Aitken (not implemented). Proper implementation requires
      % a liquid_fraction_function/meltzone_transform-consistent update of
      % T, f_ice, and f_liq.
      % T = min(Tf, alpha * T + (1 - alpha) * T_iter);
      % if use_aitken
      %    T_0 = T;
      %    for mm = 1:numel(T)
      %       T(mm) = icemodel.numerics.aitkenscalar(T_2(mm), T_1(mm), ...
      %          T_0(mm), T(mm), ...
      %          jumpmax);
      %    end
      %    T_2 = T_1;
      %    T_1 = T_0;
      % end
   end

   % Exhausting the loop is a failed solve even though the zero-based loop
   % index is maxiter-1. The caller then restores the accepted checkpoint,
   % shortens dt, and retains its bounded force-advance fallback.
   if ~ok
      iter = maxiter;
   end

   if ~ok && debug
      dumpIceEnbalFailure("maxiter_nonconvergence", T, T_ice_old, ...
         T_iter, f_ice, f_liq, f_wat, k_eff, Sc, dt, T_sfc, iM, iter, ...
         maxiter, aN, aP, aS, b);
   end

   % Subsurface energy balance linearization error [K].
   if nargout > 7
      err = icemodel.column.subsurface_linearization_error( ...
         T, T_ice_old, f_ice, f_liq, f_liq_old, dro_vapdT, ...
         dHdT, Sc, dz, dt, a2, Fc, Fp, a1);
   end

   % Surface energy balance linearization error [K]. Currently disabled, retain
   % for reference.
   % surface_err = icemodel.column.surface_linearization_error( ...
   %    T_sfc, T(1), Fc, Fp, a1); %#ok<NASGU>
end

function dumpIceEnbalFailure(reason, T, T_old, T_iter, f_ice, f_liq, ...
      f_wat, k_eff, Sc, dt, Ts, iM, iter, maxiter, aN, aP, aS, b)
   %DUMPICEENBALFAILURE Save enthalpy-solver failure diagnostics on demand.

   debug_file = getenv('ICEMODEL_DEBUG_ICEENBAL_FILE');
   if isempty(debug_file)
      return
   end

   debug_state = struct();
   debug_state.timestamp_utc = datetime('now', 'TimeZone', 'UTC');
   debug_state.reason = reason;
   debug_state.iter = iter;
   debug_state.maxiter = maxiter;
   debug_state.dt = dt;
   debug_state.Ts = Ts;
   debug_state.max_abs_dT = max(abs(T - T_iter));
   debug_state.iM = iM;
   debug_state.T = T;
   debug_state.T_old = T_old;
   debug_state.T_iter = T_iter;
   debug_state.f_ice = f_ice;
   debug_state.f_liq = f_liq;
   debug_state.f_wat = f_wat;
   debug_state.k_eff = k_eff;
   debug_state.Sc = Sc;
   debug_state.aN = aN;
   debug_state.aP = aP;
   debug_state.aS = aS;
   debug_state.b = b;

   save(debug_file, 'debug_state');
end
