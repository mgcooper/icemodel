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
         [dt_sum, subfail, OK] = NEWTIMESTEP(f_liq);

         % SURFACE TEMPERATURE
         ea = VAPPRESS(tair(metiter), Tf, liqflag) * rh(metiter) / 100;
         [Ts, ok] = SEBSOLVE(tair(metiter), swd(metiter), lwd(metiter), ...
            albedo(metiter), wspd(metiter), ppt(metiter), tppt(metiter), ...
            psfc(metiter), De(metiter), ea, cv_air, cv_liq, emiss, SB, Tf, ...
            chi, roL, scoef, liqflag, Ts, T, k_eff, dz, opts.seb_solver);
         Ts = MELTTEMP(Ts, Tf);
         if not(ok)
            sebfail_count = sebfail_count + 1;
         end
         xTs = Ts;

         while dt_sum + TINY < dt_FULL_STEP

            % HEAT CONDUCTION
            [T, OK, N] = SKINSOLVE(T, f_ice, f_liq, dz, delz, fn, dt, JJ, ...
               Ts, k_liq, cv_ice, cv_liq, ro_ice, Ls, Rv, Tf, 1e-2, opts.maxiter);

            % PROGRESS MESSAGE (SLOWS DOWN THE CODE A LOT)
            if debug == true
               % fprintf('iter = %d (%.2f%%), dt = %.0f, success = %s\n', ...
               %    iter, 100*iter/maxiter, dt, mat2str(all(OK)))
            end

            % ADAPTIVE TIME STEP
            if not(OK)
               [Ts, T, f_ice, f_liq, subfail, subiter, dt] ...
                  = RESETSUBSTEP(xTs, xT, xf_ice, xf_liq, dt_FULL_STEP, ...
                  subiter, maxsubiter, subfail, dt_sum);
               if subfail < maxsubiter
                  continue
               else
                  % disp('subfail == maxsubiter')
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
                  {Ts, Qm, Qe, Qh, Qc, chi, balance, dt_sum, ok, OK, N}, ...
                  {T, f_ice, f_liq});
            end
         end

         % MOVE TO THE NEXT TIMESTEP
         [metiter, subiter, dt] = NEXTSTEP(metiter, subiter, ...
            dt, dt_FULL_STEP, maxsubiter, OK && ok);

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
