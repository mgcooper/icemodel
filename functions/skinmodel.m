function [ice1, ice2] = skinmodel(opts) %#codegen

   % Load the physical constants
   [  cv_air, cv_liq, cv_ice, emiss, SB, epsilon, k_liq, Ls, ro_air, ...
      ro_ice, ro_liq, roLs, roLv, Rv, Tf] = ...
      icemodel.physicalConstant( ...
      'cv_air', 'cv_liq', 'cv_ice', 'emiss', 'SB', 'epsilon', 'k_liq', ...
      'Ls', 'ro_air', 'ro_ice', 'ro_liq', 'roLs', 'roLv', 'Rv', 'Tf');
   TINY = 1e-8;
   chi = 1.0;

   % Load the forcing data
   [tair, swd, lwd, albedo, wspd, rh, psfc, De, time] = METINIT(opts, 1);

   % Initialize the ice column
   [f_ice, f_liq, T, ~, ~, ~, ~, cp_sno, k_eff, dz, fn, delz, ~, ~, ~, ...
      JJ_therm, ~, ~, ~, scoef, ro_sno, ~, ~, xTsfc, xf_liq, roL, Qc, ~, ...
      fopts, ~, liqflag, ice1, ice2] = ICEINIT(opts, tair);
   
   % INITIALIZE TIMESTEPPING
   [metiter, subiter, maxiter, maxsubiter, dt_FULL_STEP, dt_min, dt_max, ...
      dt_new, numyears, ~, numspinup] = INITTIMESTEPS(opts, time);

   %% START ITERATIONS OVER YEARS
   for thisyear = 1:numyears

      iter = 1;
      while iter <= maxiter

         % INITIALIZE NEW TIMESTEP
         [dt_sum, subfail, dt_flag, OK] = INITSUBSTEP(f_liq);

         while dt_sum < dt_FULL_STEP

            % SURFACE ENERGY BALANCE
            [Qm, Qf, Qh, Qe, Qc, ~, balance, Tsfc] = ENBALANCE( ...
               tair(metiter), wspd(metiter), rh(metiter), swd(metiter), ...
               lwd(metiter), albedo(metiter), psfc(metiter), De(metiter), ...
               T, k_eff, Tf, dz, chi, xTsfc, cv_air, emiss, SB, roL, scoef, ...
               epsilon, fopts, liqflag, false);

            % HEAT CONDUCTION
            k_eff = GETGAMMA(T, f_ice, f_liq, ro_ice, k_liq, Ls, Rv, Tf);
            [T, OK] = SKINSOLVE(Tsfc, T, k_eff, ro_sno, cp_sno, dz, dt_new, ...
               JJ_therm, fn, delz, f_liq, f_ice, Tf, Rv, Ls);
            
            % ADAPTIVE TIME STEP
            if not(OK) && subiter < maxsubiter
               [subfail, subiter, dt_new, T, Tsfc, f_ice, f_liq] ...
                  = RESETSUBSTEP( xT, xTsfc, xf_ice, xf_liq, dt_max, ...
                  subiter, maxsubiter, subfail);
               continue
            else
               % UPDATE DENSITY, HEAT CAPACITY, AND SUBSTEP TIME
               k_eff = GETGAMMA(T, f_ice, f_liq, ro_ice, k_liq, Ls, Rv, Tf);
               [xT, xTsfc, xf_ice, xf_liq, dt_sum, dt_new, dt_flag, ...
                  liqflag, roL, ro_sno, cp_sno] = UPDATESUBSTEP(T, Tsfc, ...
                  f_ice, f_liq, ro_ice, ro_liq, ro_air, cv_ice, cv_liq, ...
                  dt_FULL_STEP, dt_sum, dt_new, roLv, roLs, dt_min, TINY);
            end
         end

         % SAVE OUTPUT IF SPINUP IS FINISHED
         if thisyear >= numspinup
            [ice1, ice2] = SAVEOUTPUT(iter, ice1, ice2, ...
               opts.vars1, opts.vars2, ...
               {Tsfc, Qm, Qf, Qe, Qh, Qc, chi, balance, dt_sum},  ...
               {T, f_ice, f_liq});
            % sector runs only saved ice2.T
         end

         % MOVE TO THE NEXT TIMESTEP
         [iter, metiter, subiter, dt_new] = NEXTSTEP(iter, metiter, subiter, ...
            dt_flag, dt_max, OK, dt_new);

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
