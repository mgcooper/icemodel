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
   [tair, swd, lwd, albedo, wspd, rh, psfc, ppt, tppt, De, scoef, time] ...
      = METINIT(opts, 1);

   ppt = 0 * ppt;

   % INITIALIZE THE ICE COLUMN
   [ice1, ice2, T, f_ice, f_liq, k_eff, fn, dz, delz, roL, liqflag, Ts, JJ] ...
      = ICEINIT(opts, tair);

   % INITIALIZE TIMESTEPPING
   [metiter, subiter, maxiter, maxsubiter, dt, dt_FULL_STEP, ...
      numyears, numspinup] = INITTIMESTEPS(opts, time);

   % INITIALIZE PAST VALUES
   [xT, xf_ice, xf_liq] = RESETSUBSTEP(T, f_ice, f_liq);

   %% START ITERATIONS OVER YEARS
   for thisyear = 1:numyears

      for iter = 1:maxiter

         % INITIALIZE NEW TIMESTEP
         [dt_sum, subfail, OK] = NEWTIMESTEP(f_liq);

         if ppt(metiter) > 0
            % dt = dt_FULL_STEP / maxsubiter;
         end

         % SURFACE TEMPERATURE
         ea = VAPPRESS(tair(metiter), Tf, liqflag) * rh(metiter) / 100;
         [Ts, ok] = SEBSOLVE(tair(metiter), swd(metiter), lwd(metiter), ...
            albedo(metiter), wspd(metiter), ppt(metiter), tppt(metiter), ...
            psfc(metiter), De(metiter), ea, cv_air, cv_liq, emiss, SB, Tf, ...
            chi, roL, scoef, liqflag, Ts, T, k_eff, dz, opts.seb_solver);
         Ts = MELTTEMP(Ts, Tf);

         while dt_sum + TINY < dt_FULL_STEP
            
%             [Ts, ok] = SEBSOLVE(tair(metiter), swd(metiter), lwd(metiter), ...
%                albedo(metiter), wspd(metiter), ppt(metiter), tppt(metiter), ...
%                psfc(metiter), De(metiter), ea, cv_air, cv_liq, emiss, SB, Tf, ...
%                chi, roL, scoef, liqflag, Ts, T, k_eff, dz, opts.seb_solver);
%             Ts = MELTTEMP(Ts, Tf);

            % HEAT CONDUCTION
            [T, OK, N] = SKINSOLVE(T, f_ice, f_liq, dz, delz, fn, dt, JJ, ...
               Ts, k_liq, cv_ice, cv_liq, ro_ice, Ls, Rv, Tf, 1e-2, opts.maxiter);

            % ADAPTIVE TIME STEP
            if not(OK)
               [T, f_ice, f_liq, subfail, subiter, dt] ...
                  = RESETSUBSTEP(xT, xf_ice, xf_liq, dt_FULL_STEP, ...
                  subiter, maxsubiter, subfail, dt_sum);
               if subfail < maxsubiter
                  continue
               end
            end

            % UPDATE DENSITY, HEAT CAPACITY, AND SUBSTEP TIME
            [xT, xf_ice, xf_liq, dt_sum, dt, liqflag, roL] ...
               = UPDATESUBSTEP(T, f_ice, f_liq, dt_FULL_STEP, dt_sum, ...
               dt, TINY, ro_ice, ro_liq, ro_air, cv_ice, cv_liq, roLv, roLs);
         end

         assertF(@() dt_sum < dt_FULL_STEP + 2 * TINY)

         % UPDATE SURFACE FLUXES
         k_eff = GETGAMMA(T, f_ice, f_liq, ro_ice, k_liq, Ls, Rv, Tf);
         [Qe, Qh, Qc, Qm, Qf, balance] = SEBFLUX(T, Ts, tair(metiter), ...
            swd(metiter), lwd(metiter), albedo(metiter), wspd(metiter), ...
            ppt(metiter), tppt(metiter), psfc(metiter), De(metiter), ea, ...
            Tf, k_eff, dz, cv_air, cv_liq, roL, emiss, SB, chi, epsilon, ...
            scoef, liqflag);

         % SAVE OUTPUT IF SPINUP IS FINISHED
         if thisyear >= numspinup
            % [ice1, ice2] = SAVEOUTPUT(iter, ice1, ice2, ...
            %    opts.vars1, opts.vars2, ...
            %    {Ts, Qm, Qe, Qh, Qc, chi, balance, dt_sum, ok, OK, N}, ...
            %    {T, f_ice, f_liq});
            
            % sector runs
            [ice1, ice2] = SAVEOUTPUT(iter, ice1, ice2, ...
               opts.vars1, opts.vars2, ...
               {Ts, Qm, Qe},  ...
               {T});
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
