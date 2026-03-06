function [Ts, T, f_ice, f_liq, n_subfail, substep, dt, ok] = CHECKSUBSTEP( ...
      Ts, T, f_ice, f_liq, xTs, xT, xf_ice, xf_liq, ro_ice, ro_liq, dt_sum, ...
      dt, dt_FULL_STEP, timestep, numsteps, substep, maxsubstep, n_subfail, ...
      debug, eps, ok)

   if ok
      if debug == true
         % Mass conservation / control volume check
         assertF(@() all(f_ice + f_liq * ro_liq / ro_ice <= 1 + eps))
      end
   else
      % Adapt the timestep (shorten dt and restart substep on failure)
      [Ts, T, f_ice, f_liq, n_subfail, substep, dt] ...
         = RESETSUBSTEP(xTs, xT, xf_ice, xf_liq, dt_FULL_STEP, ...
         substep, maxsubstep, n_subfail, dt_sum);

      if debug == true
         % Print progress message
         fprintf('timestep = %d (%.2f%%), dt = %.0f, ok = %s\n', ...
            timestep, 100*timestep/numsteps, dt, mat2str(ok))
      end

      % On maxsubstep, force advance using last accepted state
      if n_subfail == maxsubstep
         ok = true;
         fprintf('timestep = %d, n_subfail == maxsubstep\n', timestep)
      end
   end
