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
   opts = icemodel.configureRun(opts);
   opts = icemodel.prepareRunOutput(opts);

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

            % COUPLED SURFACE-SUBSURFACE ENERGY BALANCE
            [Ts, T, f_ice, f_liq, k_eff, ok_seb, ok_ieb, ok, n_iters] ...
               = SKINEBSOLVE(xTs, xT, xf_ice, xf_liq, dz, delz, fn, dt, JJ, ...
               ro_ice, k_liq, cv_ice, cv_liq, Ls, Rv, Tf, tair(metstep), ...
               swd(metstep), lwd(metstep), albedo(metstep), wspd(metstep), ...
               ppt(metstep), tppt(metstep), psfc(metstep), De(metstep), ...
               ea, cv_air, emiss, SB, chi, roL, scoef, liqflag, seb_solver, ...
               tol, maxiter, alpha, cpl_maxiter, cpl_Ts_tol, cpl_seb_tol, ...
               cpl_alpha, cpl_aitken, cpl_jumpmax);

            % CHECK SUBSTEP FAILURE (shorten dt and restart substep on failure)
            [Ts, T, f_ice, f_liq, n_subfail, substep, dt, ok] = ...
               CHECKSUBSTEP(Ts, T, f_ice, f_liq, xTs, xT, xf_ice, xf_liq, ...
               ro_ice, ro_liq, dt_sum, dt, dt_FULL_STEP, timestep, numsteps, ...
               substep, maxsubstep, n_subfail, debug, eps, ok);
            if not(ok)
               continue
            end

            % CHECKPOINT STATE AND SUBSTEP TIME
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

            if strcmp(opts.output_profile, 'minimal')

               [ice1, ice2] = SAVEOUTPUT(timestep, ice1, ice2, ...
                  opts.vars1, opts.vars2, ...
                  {Ts, Qm, Qe},  ...
                  {T});

            elseif strcmp(opts.output_profile, 'standard')

               [ice1, ice2] = SAVEOUTPUT(timestep, ice1, ice2, ...
                  opts.vars1, opts.vars2, ...
                  {Ts, Qm, Qe, Qh, Qc, chi, Qbal, dt_sum, ok_seb, ok_ieb, n_iters}, ...
                  {T, f_ice, f_liq});
            end
         end

         % MOVE TO THE NEXT TIMESTEP
         [metstep, substep, dt] = NEXTSTEP(metstep, substep, ...
            dt, dt_FULL_STEP, maxsubstep, ok, n_subfail, n_iters);

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
