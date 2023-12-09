function [ice1, ice2] = skinmodel(opts) %#codegen

   % INITIALIZE THE MODEL

   debug = false;
   assertF off

   % LOAD PHYSICAL CONSTANTS AND PARAMETERS
   [cv_air, cv_liq, cv_ice, emiss, SB, epsilon, k_liq, Ls, ro_air, ...
      ro_ice, ro_liq, roLs, roLv, Rv, Tf] = icemodel.physicalConstant( ...
      'cv_air', 'cv_liq', 'cv_ice', 'emiss', 'SB', 'epsilon', 'k_liq', ...
      'Ls', 'ro_air', 'ro_ice', 'ro_liq', 'roLs', 'roLv', 'Rv', 'Tf');
   TINY = 1e-8;
   chi = 1.0;

   % LOAD THE FORCING DATA
   [tair, swd, lwd, albedo, wspd, rh, psfc, De, scoef, time] = METINIT(opts, 1);

   % INITIALIZE THE ICE COLUMN
   [ice1, ice2, T, f_ice, f_liq, k_eff, fn, dz, delz, roL, liqflag, Ts, JJ] ...
      = ICEINIT(opts, tair);

   % INITIALIZE TIMESTEPPING
   [metiter, subiter, maxiter, maxsubiter, dt, dt_FULL_STEP, ...
      numyears, numspinup] = INITTIMESTEPS(opts, time);

   %% START ITERATIONS OVER YEARS
   for thisyear = 1:numyears

      for iter = 1:maxiter

         % INITIALIZE NEW TIMESTEP
         [dt_sum, subfail, OK] = NEWTIMESTEP(f_liq);

         while dt_sum + TINY < dt_FULL_STEP

            % SURFACE ENERGY BALANCE
            [Qm, Qf, Qh, Qe, Qc, ~, balance, Ts] = ENBALANCE( ...
               tair(metiter), swd(metiter), lwd(metiter), albedo(metiter), ...
               wspd(metiter), rh(metiter), psfc(metiter), De(metiter), ...
               T, k_eff, Tf, dz, chi, xTs, cv_air, emiss, SB, roL, scoef, ...
               epsilon, liqflag, false);

            % HEAT CONDUCTION
            k_eff = GETGAMMA(T, f_ice, f_liq, ro_ice, k_liq, Ls, Rv, Tf);
            [T, OK] = SKINSOLVE(T, f_ice, f_liq, ro_sno, cp_sno, k_eff, fn, ...
               delz, dz, dt, JJ, Ts, Tf, Rv, Ls);

            % ADAPTIVE TIME STEP
            if not(OK)
               [T, Ts, f_ice, f_liq, subfail, subiter, dt] ...
                  = RESETSUBSTEP(xT, xTs, xf_ice, xf_liq, dt_FULL_STEP, ...
                  subiter, maxsubiter, subfail, dt_sum);
               if subfail < maxsubiter
                  continue
               end
            end

            % UPDATE DENSITY, HEAT CAPACITY, AND SUBSTEP TIME
            [xT, xTs, xf_ice, xf_liq, dt_sum, dt, liqflag, roL] ...
               = UPDATESUBSTEP(T, Ts, f_ice, f_liq, dt_FULL_STEP, dt_sum, ...
               dt, TINY, ro_ice, ro_liq, ro_air, cv_ice, cv_liq, roLv, roLs);
         end

         assertF(@() dt_sum < dt_FULL_STEP + 2 * TINY)

         % SAVE OUTPUT IF SPINUP IS FINISHED
         if thisyear >= numspinup
            [ice1, ice2] = SAVEOUTPUT(iter, ice1, ice2, ...
               opts.vars1, opts.vars2, ...
               {Ts, Qm, Qf, Qe, Qh, Qc, chi, balance, dt_sum},  ...
               {T, f_ice, f_liq});
            % sector runs only saved ice2.T
         end

         % MOVE TO THE NEXT TIMESTEP
         [metiter, subiter, dt] = NEXTSTEP(metiter, subiter, ...
            dt, dt_FULL_STEP, maxsubiter, OK);

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
end
