function [ice1, ice2] = skinmodel(opts)
   %SKINMODEL Surface energy balance model for glacier ice.
   %
   %
   % See also: icemodel
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
   sebfail_count = 0;
   maxcpliter = opts.maxcpliter;
   cpltol = opts.cpltol;
   omega = opts.omega;
   sebtol = opts.sebtol;
   use_aitken = opts.use_aitken;
   aitken_jumpmax = opts.aitken_jumpmax;

   % LOAD THE FORCING DATA
   [tair, swd, lwd, albedo, wspd, rh, psfc, ppt, tppt, De, scoef, time] ...
      = METINIT(opts, 1);

   % INITIALIZE THE ICE COLUMN
   [ice1, ice2, T, f_ice, f_liq, k_eff, fn, dz, delz, roL, liqflag, Ts, JJ] ...
      = ICEINIT(opts, tair);

   % INITIALIZE TIMESTEPPING
   [metiter, subiter, maxiter, maxsubiter, dt, dt_FULL_STEP, ...
      numyears, numspinup] = INITTIMESTEPS(opts, time);

   % INITIALIZE PAST VALUES
   [xTs, xT, xf_ice, xf_liq] = RESETSUBSTEP(Ts, T, f_ice, f_liq);

   %% START ITERATIONS OVER YEARS
   for thisyear = 1:numyears

      for iter = 1:maxiter

         % INITIALIZE NEW TIMESTEP
         [dt_sum, subfail, ok_seb, ok_ieb] = NEWTIMESTEP(f_liq, opts);

         % Atmospheric vapor pressure (fixed over this full forcing step).
         ea = VAPPRESS(tair(metiter), Tf, liqflag) * rh(metiter) / 100;

         while dt_sum + TINY < dt_FULL_STEP

            % Pre-coupler Ts predictor using checkpoint state
            k_eff = GETGAMMA(xT, xf_ice, xf_liq, ro_ice, k_liq, Ls, Rv, Tf);
            [Ts, ok_seb] = SEBSOLVE(tair(metiter), swd(metiter), ...
               lwd(metiter), albedo(metiter), wspd(metiter), ...
               ppt(metiter), tppt(metiter), psfc(metiter), De(metiter), ...
               ea, cv_air, cv_liq, emiss, SB, Tf, chi, roL, scoef, ...
               liqflag, xTs, xT, k_eff, dz, opts.seb_solver);
            Ts = MELTTEMP(Ts, Tf);

            % Initial values for Aitken acceleration
            Ts_1 = nan;
            Ts_2 = nan;

            % Run the coupler
            ok_cpl = false;
            for cpliter = 1:maxcpliter

               % Solve subsurface from checkpoint state w/o physical advancement
               [T, ok_ieb, N] = SKINSOLVE(xT, xf_ice, xf_liq, dz, delz, fn, ...
                  dt, JJ, Ts, k_liq, cv_ice, cv_liq, ro_ice, Ls, Rv, Tf, ...
                  1e-2, opts.maxiter);
               if not(ok_ieb)
                  break
               end

               % Solve surface (in-loop corrector using updated trial state)
               old = Ts;
               k_eff = GETGAMMA(T, f_ice, f_liq, ro_ice, k_liq, Ls, Rv, Tf);
               [Ts, ok_seb] = SEBSOLVE(tair(metiter), swd(metiter), ...
                  lwd(metiter), albedo(metiter), wspd(metiter), ...
                  ppt(metiter), tppt(metiter), psfc(metiter), De(metiter), ...
                  ea, cv_air, cv_liq, emiss, SB, Tf, chi, roL, scoef, ...
                  liqflag, Ts, T, k_eff, dz, opts.seb_solver);
               Ts = MELTTEMP(Ts, Tf);

               if not(ok_seb)
                  sebfail_count = sebfail_count + 1;
                  break
               end

               % SEB residual
               sebres = abs(fSEB(Ts, tair(metiter), swd(metiter), ...
                  lwd(metiter), albedo(metiter), wspd(metiter), ...
                  ppt(metiter), tppt(metiter), psfc(metiter), ...
                  De(metiter), ea, cv_air, cv_liq, emiss, SB, Tf, ...
                  chi, roL, scoef, CONDUCT(k_eff, T, dz, Ts), liqflag));

               % At melt cap (Ts ~= Tf), positive residual is melt energy (Qm).
               if Ts >= Tf
                  sebres = 0.0;
               end

               % Check convergence
               if abs(Ts - old) < cpltol && sebres < sebtol
                  ok_cpl = true;
                  break
               end

               % Aitken acceleration with relaxation-fallback
               Ts_0 = Ts;
               Ts = MELTTEMP(aitkenscalar(Ts_2, Ts_1, Ts_0, ...
                  (1.0 - omega) * old + omega * Ts, ... % relaxation fallback
                  aitken_jumpmax, use_aitken), Tf);
               Ts_2 = Ts_1;
               Ts_1 = Ts_0;
            end
            OK = ok_seb && ok_ieb && ok_cpl;

            % PROGRESS MESSAGE (SLOWS DOWN THE CODE A LOT)
            if debug == true && not(OK)
               fprintf('iter = %d (%.2f%%), dt = %.0f, success = %s\n', ...
                  iter, 100*iter/maxiter, dt, mat2str(all(ok_ieb)))
            end

            % ADAPTIVE TIME STEP
            if not(OK)
               [Ts, T, f_ice, f_liq, subfail, subiter, dt] ...
                  = RESETSUBSTEP(xTs, xT, xf_ice, xf_liq, dt_FULL_STEP, ...
                  subiter, maxsubiter, subfail, dt_sum);
               if subfail < maxsubiter
                  continue
               else
                  fprintf('iter = %d, subfail == maxsubiter\n', iter)
               end
            end

            % UPDATE DENSITY, HEAT CAPACITY, AND SUBSTEP TIME
            [xTs, xT, xf_ice, xf_liq, dt_sum, dt, liqflag, roL] ...
               = UPDATESUBSTEP(Ts, T, f_ice, f_liq, dt_FULL_STEP, dt_sum, ...
               dt, TINY, ro_ice, ro_liq, ro_air, cv_ice, cv_liq, roLv, roLs);
         end

         assertF(@() dt_sum < dt_FULL_STEP + 2 * TINY)

         % UPDATE SURFACE FLUXES
         k_eff = GETGAMMA(T, f_ice, f_liq, ro_ice, k_liq, Ls, Rv, Tf);
         [Qe, Qh, Qc, Qm, ~, balance] = SEBFLUX(T, Ts, tair(metiter), ...
            swd(metiter), lwd(metiter), albedo(metiter), wspd(metiter), ...
            ppt(metiter), tppt(metiter), psfc(metiter), De(metiter), ea, ...
            Tf, k_eff, dz, cv_air, cv_liq, roL, emiss, SB, chi, epsilon, ...
            scoef, liqflag);

         % SAVE OUTPUT IF SPINUP IS FINISHED
         if thisyear >= numspinup

            if strcmp('sector', opts.sitename)
               [ice1, ice2] = SAVEOUTPUT(iter, ice1, ice2, ...
                  opts.vars1, opts.vars2, ...
                  {Ts, Qm, Qe},  ...
                  {T});
            else
               [ice1, ice2] = SAVEOUTPUT(iter, ice1, ice2, ...
                  opts.vars1, opts.vars2, ...
                  {Ts, Qm, Qe, Qh, Qc, chi, balance, dt_sum, ok_seb, ok_ieb, N}, ...
                  {T, f_ice, f_liq});
            end
         end

         % MOVE TO THE NEXT TIMESTEP
         [metiter, subiter, dt] = NEXTSTEP(metiter, subiter, ...
            dt, dt_FULL_STEP, maxsubiter, ok_ieb, subfail, N);

      end % timesteps (one year)

      % RESTART THE MET DATA ITERATOR DURING SPIN UP
      if thisyear < numspinup
         metiter = 1;
         continue
      end

      % WRITE TO DISK
      WRITEOUTPUT(ice1, ice2, opts, thisyear-numspinup+1, ...
         time((thisyear-numspinup)*maxiter+1:(thisyear-numspinup+1)*maxiter), ...
         swd, lwd, albedo)
   end
   fprintf('SEBSOLVE fail count=%d\n', sebfail_count)
end
