function [T, f_ice, f_liq, ok] = MZTRANSFORM(T, T_old, f_liq, f_wat, ...
      ro_ice, ro_liq, Tf, TL, TH, fcp, f_liq_min, f_liq_max, i_M, ok)
   %MZTRANSFORM Liquid fraction change to temperature-enthalpy transformation
   %
   % This function uses the change in liquid fraction returned by the numerical
   % solution of the enthalpy equation to update the temperature of nodes
   % undergoing phase change, using the temperature-enthalpy relationship.
   % This ensures the change in enthalpy due to solid-liquid phase change is
   % accounted for in terms of both the latent heat and the specific heat of the
   % control volume, i.e., that the temperature, liquid fraction, and enthalpy
   % are consistent. This is sometimes referred to as a "corrector" step. Note
   % that for nodes undergoing phase change, the input to this function, T, is
   % the change in liquid fraction due to solid-liquid phase change multiplied
   % by the density of liquid water. For nodes which are not undergoing phase
   % change, T is the temperature.
   %
   % This function also implements three numerical checks:
   %  1) if the enthalpy predictor implies an impossible liquid fraction
   %     (negative, above the available water, or above the available ice),
   %  2) if a melt-zone node overshoots either melt-zone phase boundary by
   %     more than 5%, meaning the local melt-curve linearization is no longer
   %     trusted, and
   %  3) if a node that started outside the melt zone skipped it entirely,
   %     meaning a node below TL ended the step above TH (or vice versa).
   %
   % If any are true, OK is set false, control returns to the main program,
   % the time step is shortened, and the solution is repeated until a valid
   % solution is obtained.
   %
   % Nodes that already started inside the melt zone are allowed to leave it
   % during the corrector step as long as the predictor overshoots the melt-
   % zone phase bounds only slightly. Small overshoots are projected onto the
   % correct branch using the old melt-zone slope; larger overshoots trigger
   % timestep reduction.
   %
   % Notes:
   %  - T_old is needed both for the melt-zone skipping check and for the
   %    linearized projection of nodes that exit the melt zone.
   %  - f_liq_min/max are f_liq at T=TL/TH, given the current f_wat, not the
   %    min/max possible f_liq (but f_liq_max = f_wat for an ice model)
   %  - i_M is the melt zone nodes at each inner iteration, f_wat is the
   %    water fraction at the start of the outer iterations.
   %
   % See also: FREEZECURVE MELTCURVE
   %
   %#codegen

   % Early exit if the predictor implies a liquid-fraction increase larger
   % than the maximum liquid-fraction-equivalent ice.
   if any((T(i_M) / ro_liq) > (ro_ice / ro_liq - f_liq(i_M)))
      ok = false;
      return
   end

   % Overshoot tolerance
   tol = 0.05;
   sqrt_eps = sqrt(eps);

   %%% Update liquid / solid fractions
   %
   % Note that here, T(i_M) = ro_liq * (f_liq_new - f_liq_old), T(~i_M) = T_new

   % Update the liquid fraction of melt-zone layers (f_liq = f_liq_old + P/ro)
   f_liq(i_M) = f_liq(i_M) + T(i_M) / ro_liq; % line 79 of ftemp.f

   % This is what the transformation above does:
   %
   % f_liq_new = P / ro_liq + f_liq_old
   %           = ro_liq * (f_liq_new - f_liq_old) / ro_liq + f_liq_old
   %           =          (f_liq_new - f_liq_old)          + f_liq_old
   % f_liq_new = f_liq_new

   % Update the ice fraction implied by the current predictor state.
   f_ice = (f_wat - f_liq) * ro_liq / ro_ice;

   % Jordan shortens the timestep if Pmelt overshoots either melt-zone
   % boundary by more than 5%. Express the same check in the current
   % volumetric-liquid-fraction variables.
   if any( ...
         (f_liq(i_M) < (1.0 - tol) * f_liq_min(i_M)) | ...
         (f_liq(i_M) > (1.0 + tol) * f_liq_max(i_M)) ...
         )

      dumpMZTransformFailure("pmelt_overshot_melt_zone_bounds", T, T_old, ...
         f_ice, f_liq, f_wat, f_liq_min, f_liq_max, i_M);

      ok = false;
      return
   end

   %%% Update temperature
   % Below here, transform T(i_M) = ro_liq * (f_liq_new - f_liq_old) -> T_new
   %
   % The older version applied the melt-zone analytic inverse to all melt-zone
   % nodes that survived the tol overshoot check (like Jordan). That is not
   % correct for nodes that exited the melt zone during the predictor step.
   % Below treats nodes that remain in the melt zone separately from nodes that
   % exit it.

   % Identify nodes that left the melt zone during the predictor step but not
   % by enough to reject the timestep. These transitions are handled below by
   % projecting the predictor through the local melt-curve slope and then
   % letting the next Picard iteration solve on the correct melt curve branch.
   i_lo = false(size(i_M));
   i_hi = false(size(i_M));
   i_lo(i_M) = f_liq(i_M) < f_liq_min(i_M);
   i_hi(i_M) = f_liq(i_M) > f_liq_max(i_M);

   % A melt-zone node may legitimately leave the phase-change interval during
   % the corrector step. For nodes that remain within the melt zone, use the
   % analytic inverse phase function with the new f_liq to update T:
   i_ok = i_M & ~(i_lo | i_hi);
   if any(i_ok)
      T(i_ok) = Tf - sqrt(f_wat(i_ok) ./ f_liq(i_ok) - 1.0) / fcp; % Eq. 133a
   end

   % Note: above uses the new f_liq directly, which is correct. Don't use this
   % update (use it below for nodes that have exited):
   % T(i_M) = T_old(i_M) + T(i_M) ./ (ro_liq * dFdT_old(i_M));

   % For nodes that exit the melt zone with only a small predictor overshoot,
   % use the local linearized predictor dFdT to move T onto the frozen or liquid
   % branch, then let the next nonlinear iteration rebuild the coefficients with
   % the updated classification.

   % At this point, T(i_ok) has already been updated to temperature, but
   % T(i_lo | i_hi) still stores P = ro_liq * (f_liq_new - f_liq_old).
   % The projection below must therefore only be used for the exit nodes.
   if any(i_lo | i_hi)

      % What below does conceptually:
      % T_new ≈ T_old + P / (ro_liq * (dF/dT)_old)
      %       = T_old + (f_liq_new - f_liq_old) / (dF/dT)_old
      %
      % This uses the old melt-zone slope because the predictor was formed
      % using the old melt-zone linearization, and nodes that have just exited
      % the melt zone do not yet have a branch-consistent new temperature from
      % which to evaluate an updated slope.

      % Evaluate dFdT at nodes within the melt zone using T_old
      dFdT_old = FREEZECURVE(T_old(i_M), ro_ice, ro_liq, fcp, Tf, ...
         [], [], f_wat(i_M));

      % Linearized temperature update evaluated on melt-zone nodes using
      % T_old and dFdT_old
      T_new = T_old(i_M) + T(i_M) ./ (ro_liq * dFdT_old);

      % Update nodes that exited the melt zone onto the frozen branch
      if any(i_lo(i_M))
         T(i_lo) = min(T_new(i_lo(i_M)), TL - sqrt_eps);
      end

      % Update nodes that exited the melt zone onto the liquid branch
      if any(i_hi(i_M))
         T(i_hi) = max(T_new(i_hi(i_M)), TH + sqrt_eps);
      end
   end

   % Fully crossing the melt zone in one step is a failure
   if any( ...
         T_old(~i_M) < TL & T(~i_M) >= TH) || ...
         any( ...
         T_old(~i_M) > TH & T(~i_M) < TL) || ...
         any( ...
         iscomplex(T)) || any(~isfinite(T) ...
         )

      dumpMZTransformFailure("skipped_melt_zone_or_complex", T, T_old, ...
         f_ice, f_liq, f_wat, f_liq_min, f_liq_max, i_M);

      ok  = false;
      return
   end

   % Project T/f_ice/f_liq onto the enthalpy-temperature curve (corrector step)
   [T, f_ice, f_liq] = MELTCURVE(T, f_ice, f_liq, ro_ice, ro_liq, fcp, Tf);
end

function dumpMZTransformFailure(reason, T, T_old, f_ice, f_liq, ...
      f_wat, f_liq_min, f_liq_max, i_M)
   %DUMPMZTRANSFORMFAILURE Save melt-zone failure diagnostics on demand.

   debug_file = getenv('ICEMODEL_DEBUG_MZTRANSFORM_FILE');
   if isempty(debug_file)
      return
   end

   debug_state = struct();
   debug_state.timestamp_utc = datetime('now', 'TimeZone', 'UTC');
   debug_state.reason = reason;
   debug_state.T = T;
   debug_state.T_old = T_old;
   debug_state.f_ice = f_ice;
   debug_state.f_liq = f_liq;
   debug_state.f_wat = f_wat;
   debug_state.f_liq_min = f_liq_min;
   debug_state.f_liq_max = f_liq_max;
   debug_state.i_M = i_M;

   save(debug_file, 'debug_state');
end
