function [ice1, ice2] = skinmodel(opts)
   %SKINMODEL Surface energy balance model for glacier ice.
   %
   %
   % See also: icemodel, icemodel.setopts
   %
   %#codegen

   %% INITIALIZE THE MODEL

   debug = true;
   assertF off

   % LOAD PHYSICAL CONSTANTS AND PARAMETERS
   [cv_air, cv_ice, cv_liq, emiss, SB, epsilon, k_liq, Ls, ro_air, ...
      ro_ice, ro_liq, roLs, roLv, Rv, Tf] = icemodel.physicalConstant( ...
      'cv_air', 'cv_ice', 'cv_liq', 'emiss', 'SB', 'epsilon', 'k_liq', ...
      'Ls', 'ro_air', 'ro_ice', 'ro_liq', 'roLs', 'roLv', 'Rv', 'Tf');
   TINY = 1e-8;
   chi = 1.0;

   % LOAD THE FORCING DATA
   [tair, swd, lwd, albedo, wspd, rh, psfc, ppt, tppt, De, scoef, time] ...
      = METINIT(opts, 1);

   % INITIALIZE THE THERMAL MODEL
   [ice1, ice2, T, f_ice, f_liq, k_eff, fn, dz, delz, roL, liqflag, Ts, JJ] ...
      = ICEINIT(opts, tair);

   % INITIALIZE TIMESTEPPING
   [metstep, substep, numsteps, maxsubstep, dt, dt_FULL_STEP, ...
      numyears, numspinup] = INITTIMESTEPS(opts, time);

   % INITIALIZE PAST VALUES
   [xTs, xT, xf_ice, xf_liq] = RESETSUBSTEP(Ts, T, f_ice, f_liq);

   % UNPACK SOLVER OPTS
   [bc_type, seb_solver, cpl_maxiter, cpl_Ts_tol, cpl_seb_tol, ...
      cpl_alpha, cpl_aitken, cpl_jumpmax, maxiter, tol, alpha, ...
      use_aitken, jumpmax] = icemodel.getopts(opts, ...
      'bc_type', 'seb_solver', 'cpl_maxiter', 'cpl_Ts_tol', ...
      'cpl_seb_tol', 'cpl_alpha', 'cpl_aitken', 'cpl_jumpmax', ...
      'maxiter', 'tol', 'alpha', 'use_aitken', 'jumpmax');

   %% START TIMESTEPS OVER YEARS
   for thisyear = 1:numyears

      for timestep = 1:numsteps

         % INITIALIZE NEW TIMESTEP
         [dt_sum, n_subfail, ok_seb, ok_ieb] = NEWTIMESTEP(f_liq, bc_type);

         % SURFACE TERMS (atmospheric vapor pressure fixed over this full step)
         ea = VAPPRESS(tair(metstep), Tf, liqflag) * rh(metstep) / 100;

         while dt_sum + TINY < dt_FULL_STEP

            % Pre-coupler Ts predictor using checkpoint state
            k_eff = GETGAMMA(xT, xf_ice, xf_liq, ro_ice, k_liq, Ls, Rv, Tf);
            [Ts, ok_seb] = SEBSOLVE(tair(metstep), swd(metstep), ...
               lwd(metstep), albedo(metstep), wspd(metstep), ...
               ppt(metstep), tppt(metstep), psfc(metstep), De(metstep), ...
               ea, cv_air, cv_liq, emiss, SB, Tf, chi, roL, scoef, ...
               liqflag, xTs, xT, k_eff, dz, seb_solver);
            Ts = MELTTEMP(Ts, Tf);

            % Initial values for Aitken acceleration
            Ts_1 = nan;
            Ts_2 = nan;

            % Run the coupler
            ok_cpl = false;
            for cpliter = 1:cpl_maxiter

               % Solve subsurface from checkpoint state w/o physical advancement
               [T, f_ice, f_liq, k_eff, ok_ieb, n_iters] = SKINSOLVE(xT, ...
                  xf_ice, xf_liq, dz, delz, fn, dt, JJ, Ts, k_liq, cv_ice, ...
                  cv_liq, ro_ice, Ls, Rv, Tf, tol, maxiter, alpha);
               if not(ok_ieb)
                  break
               end

               % Solve surface (in-loop corrector using updated trial state)
               old = Ts;
               [Ts, ok_seb] = SEBSOLVE(tair(metstep), swd(metstep), ...
                  lwd(metstep), albedo(metstep), wspd(metstep), ...
                  ppt(metstep), tppt(metstep), psfc(metstep), De(metstep), ...
                  ea, cv_air, cv_liq, emiss, SB, Tf, chi, roL, scoef, ...
                  liqflag, Ts, T, k_eff, dz, seb_solver);
               Ts = MELTTEMP(Ts, Tf);

               if not(ok_seb)
                  break
               end

               % SEB residual
               seb_res = abs(fSEB(Ts, tair(metstep), swd(metstep), ...
                  lwd(metstep), albedo(metstep), wspd(metstep), ...
                  ppt(metstep), tppt(metstep), psfc(metstep), ...
                  De(metstep), ea, cv_air, cv_liq, emiss, SB, Tf, ...
                  chi, roL, scoef, CONDUCT(k_eff, T, dz, Ts), liqflag));

               % At melt cap (Ts ~= Tf), positive residual is melt energy (Qm).
               if Ts >= Tf
                  seb_res = 0.0;
               end

               % Check convergence
               if abs(Ts - old) < cpl_Ts_tol && seb_res < cpl_seb_tol
                  ok_cpl = true;
                  break
               end

               % Aitken acceleration with relaxation-fallback
               Ts_0 = Ts;
               Ts = MELTTEMP(aitkenscalar(Ts_2, Ts_1, Ts_0, ...
                  (1.0 - cpl_alpha) * old + cpl_alpha * Ts, ... % relaxation
                  cpl_jumpmax, cpl_aitken), Tf);
               Ts_2 = Ts_1;
               Ts_1 = Ts_0;
            end
            ok = ok_seb && ok_ieb && ok_cpl;

            % PROGRESS MESSAGE (slows down the code a lot)
            if debug == true && not(ok)
               fprintf('timestep = %d (%.2f%%), dt = %.0f, success = %s\n', ...
                  timestep, 100*timestep/numsteps, dt, mat2str(ok))
            end

            % ADAPTIVE TIME STEP (shorten dt and restart substep on failure)
            if not(ok)
               [Ts, T, f_ice, f_liq, n_subfail, substep, dt] ...
                  = RESETSUBSTEP(xTs, xT, xf_ice, xf_liq, dt_FULL_STEP, ...
                  substep, maxsubstep, n_subfail, dt_sum);
               if n_subfail < maxsubstep
                  continue
               else
                  fprintf('timestep = %d, n_subfail == maxsubstep\n', timestep)
               end
            end

            % UPDATE STATE AND SUBSTEP TIME (checkpoint state)
            [xTs, xT, xf_ice, xf_liq, dt_sum, dt, liqflag, roL] ...
               = UPDATESUBSTEP(Ts, T, f_ice, f_liq, dt_FULL_STEP, dt_sum, ...
               dt, TINY, ro_ice, ro_liq, ro_air, cv_ice, cv_liq, roLv, roLs);
         end

         % Error if dt accumulation exceeds full step
         assertF(@() dt_sum < dt_FULL_STEP + 2 * TINY)

         % DIAGNOSE SURFACE FLUXES
         [Qe, Qh, Qc, Qm, ~, Qbal] = SEBFLUX(T, Ts, tair(metstep), ...
            swd(metstep), lwd(metstep), albedo(metstep), wspd(metstep), ...
            ppt(metstep), tppt(metstep), psfc(metstep), De(metstep), ea, ...
            Tf, k_eff, dz, cv_air, cv_liq, roL, emiss, SB, chi, epsilon, ...
            scoef, liqflag);

         % SAVE OUTPUT IF SPINUP IS FINISHED
         if thisyear >= numspinup

            if strcmp('sector', opts.sitename)

               [ice1, ice2] = SAVEOUTPUT(timestep, ice1, ice2, ...
                  opts.vars1, opts.vars2, ...
                  {Ts, Qm, Qe},  ...
                  {T});
            else
               [ice1, ice2] = SAVEOUTPUT(timestep, ice1, ice2, ...
                  opts.vars1, opts.vars2, ...
                  {Ts, Qm, Qe, Qh, Qc, chi, Qbal, dt_sum, ok_seb, ok_ieb, n_iters}, ...
                  {T, f_ice, f_liq});
            end
         end

         % MOVE TO THE NEXT TIMESTEP
         [metstep, substep, dt] = NEXTSTEP(metstep, substep, ...
            dt, dt_FULL_STEP, maxsubstep, ok_ieb, n_subfail, n_iters);

      end % timesteps (one year)

      % RESTART THE MET DATA STEP INDEX DURING SPIN UP
      if thisyear < numspinup
         metstep = 1;
         continue
      end

      % WRITE TO DISK
      WRITEOUTPUT(ice1, ice2, opts, thisyear-numspinup+1, ...
         time((thisyear-numspinup)*numsteps+1:(thisyear-numspinup+1)*numsteps), ...
         swd, lwd, albedo)
   end
end
