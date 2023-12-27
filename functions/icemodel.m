function [ice1, ice2] = icemodel(opts) %#codegen
   % ICEMODEL Simulate the phase change process in glacier ice.
   %
   % This function models the phase change process in melting glacier ice. It
   % uses iterative processes to update the temperature, liquid and ice fraction
   % at each time step. This model considers both the surface and subsurface
   % energy balance to simulate the phase change.
   %
   % Syntax:
   % [ice1, ice2] = ICEMODEL(opts)
   %
   % Inputs:
   % opts - A structure containing model options and parameters. Defined by the
   %        icemodel.setopts function. The fields include:
   %        * sitename         - Site name for the model simulation
   %        * simmodel         - Simulation model identifier
   %        * simyears         - Years for which the simulation is done
   %        * forcings         - Type of forcing data used
   %        * userdata         - User-defined data type
   %        * uservars         - User-defined variables
   %        * savedata         - Flag indicating if data should be saved
   %        * testname         - Name of the test (default: 'none')
   %        * testpoint        - Test point if applicable (default: 'none')
   %        * ... (other parameters related to the model configuration)
   %
   % Outputs:
   % ice1  - 1-dimensional data storing variables defined at the ice surface or
   %         near-surface atmosphere. Contains one value per timestep.
   % ice2  - 2-dimensional data storing variables defined on the subsurface ice
   %         column control volume mesh. Contains one column per timestep.
   % met   - Struct containing the meteorological data used for the simulation.
   %
   % See also: ICEMODEL.SETOPTS

   %% INITIALIZE THE MODEL

   debug = false;
   assertF on

   % LOAD PHYSICAL CONSTANTS AND PARAMETERS
   [cv_air, cv_liq, cv_ice, emiss, SB, epsilon, fcp, k_liq, Lf, ...
      Ls, Lv, ro_air, ro_ice, ro_liq, roLf, roLs, roLv, Rv, Tf] ...
      = icemodel.physicalConstant( ...
      'cv_air','cv_liq','cv_ice','emiss', 'SB','epsilon','fcp','k_liq','Lf', ...
      'Ls','Lv','ro_air','ro_ice','ro_liq','roLf','roLs','roLv', 'Rv','Tf');
   TINY = 1e-8;

   % LOAD THE FORCING DATA
   [tair, swd, lwd, albedo, wspd, rh, psfc, De, scoef, time] = METINIT(opts, 1);

   % INITIALIZE THE ICE COLUMN
   [ice1, ice2, T, f_ice, f_liq, k_eff, fn, dz, delz, roL, liqflag, Ts, ...
      JJ, z_therm, dz_therm, ~, Sp, Fc, Fp, TL, TH, flmin, flmax, f_min, ...
      liqresid, ro_iwe, ro_wie] = ICEINIT(opts, tair);

   % INITIALIZE THE SPECTRAL EXTINCTION COEFFICIENTS
   [I0, z_spect, spect_N, spect_S, solardwavl, dz_spect] ...
      = EXTCOEFSINIT(opts, ro_ice);

   % INITIALIZE TIMESTEPPING
   [metiter, subiter, maxiter, maxsubiter, dt, dt_FULL_STEP, ...
      numyears, numspinup] = INITTIMESTEPS(opts, time);

   % INITIALIZE PAST VALUES
   [xT, xTs, xf_ice, xf_liq] = RESETSUBSTEP(T, Ts, f_ice, f_liq);
   
   %% START ITERATIONS OVER YEARS
   for thisyear = 1:numyears

      for iter = 1:maxiter

         % INITIALIZE NEW TIMESTEP
         [dt_sum, subfail, OK, d_liq, d_drn, d_evp] = NEWTIMESTEP(f_liq);

         % update the upper layer diffusion length scale
         % if sqrt(k_eff(1)*dt/(ro_sno(1)*cp_sno(1))) > dz(1)
         %    pause;
         % end

         % SUBSURFACE SOLAR RADIATION SOURCE-TERM
         [Sc, chi] = UPDATEEXTCOEFS(swd(metiter), albedo(metiter), I0, ...
            dz, z_spect, dz_spect, z_therm, dz_therm, spect_N, spect_S, ...
            solardwavl, ro_liq * f_liq + ro_ice * f_ice);
         % icemodel.verifyShortwaveBalance(Qsi, Sc, albedo, chi, dz, 1);

         % SURFACE TEMPERATURE
         ea = VAPPRESS(tair(metiter), Tf, liqflag) * rh(metiter) / 100;
         Ts = SEBSOLVE(tair(metiter), swd(metiter), lwd(metiter), ...
            albedo(metiter), wspd(metiter), psfc(metiter), De(metiter), ...
            ea, cv_air, emiss, SB, Tf, chi, roL, scoef, liqflag, ...
            Ts, T, k_eff, dz, opts.sebsolver);

         while dt_sum + TINY < dt_FULL_STEP

            % SUBSURFACE ENERGY BALANCE
            [T, f_ice, f_liq, OK, errT, errH] = ICEENBAL(T, f_ice, f_liq, ...
               dz, delz, fn, Sc, dt, JJ, Ts, k_liq, cv_ice, cv_liq, ...
               ro_ice, ro_liq, ro_air, Ls, Lf, roLf, Rv, Tf, fcp, TL, TH, ...
               flmin, flmax, ro_iwe, ro_wie);

            if debug == true
               % ERROR MESSAGE (SLOWS DOWN THE CODE A LOT)
               fprintf('iter = %d (%.2f%%), dt = %.0f, success = %s\n', ...
                  iter,100*iter/maxiter,dt,mat2str(all(OK)))
               assertF(@() all(f_ice + f_liq * ro_wie <= 1 + eps))
            end

            % PHASE BOUNDARY OVERSHOOT, DECREASE THE TIME STEP AND START OVER
            if not(OK)
               [T, Ts, f_ice, f_liq, subfail, subiter, dt] ...
                  = RESETSUBSTEP(xT, xTs, xf_ice, xf_liq, dt_FULL_STEP, ...
                  subiter, maxsubiter, subfail, dt_sum);
               if subfail < maxsubiter
                  continue
               end
            end

            % UPDATE SURFACE FLUXES
            k_eff = GETGAMMA(T, f_ice, f_liq, ro_ice, k_liq, Ls, Rv, Tf);
            [Qe, Qh, Qc, Qm, Qf, balance] = SEBFLUX(T, Ts, tair(metiter), ...
               swd(metiter), lwd(metiter), albedo(metiter), wspd(metiter), ...
               psfc(metiter), De(metiter), ea, Tf, k_eff, dz, cv_air, roL, ...
               emiss, SB, chi, epsilon, scoef, liqflag);

            % COMPUTE MASS BALANCE
            [T, f_ice, f_liq, d_liq, d_evp, d_drn] = ICEMF(T, f_ice, f_liq, ...
               ro_ice, ro_liq, cv_ice, cv_liq, Lf, Ls, Lv, Tf, TL, fcp, ...
               xf_liq, Sc, Sp, JJ, f_min, dz_therm, dt, Qe, ro_iwe, ...
               d_liq, d_drn, d_evp, flmin, liqresid);

            % UPDATE DENSITY, HEAT CAPACITY, AND SUBSTEP TIME
            [xT, xTs, xf_ice, xf_liq, dt_sum, dt, liqflag, roL] ...
               = UPDATESUBSTEP(T, Ts, f_ice, f_liq, dt_FULL_STEP, dt_sum, ...
               dt, TINY, ro_ice, ro_liq, ro_air, cv_ice, cv_liq, roLv, roLs);

            % zD = sqrt(k_eff(1)*dt/(ro_sno(1)*cp_sno(1)));
            % if zD > dz(1) 
            % end
         end

         assertF(@() dt_sum < dt_FULL_STEP + 2 * TINY)

         % SAVE OUTPUT IF SPINUP IS FINISHED
         if thisyear >= numspinup
            [ice1, ice2] = SAVEOUTPUT(iter, ice1, ice2, ...
               opts.vars1, opts.vars2, ...
               {Ts, Qm, Qf, Qe, Qh, Qc, chi, balance, dt_sum, ok, OK, N}, ...
               {T, f_ice, f_liq, d_liq, d_drn, d_evp, Sc});

            % For a regional run:
            % {Ts},  ...
            % {T, f_ice, f_liq, d_liq, d_drn, d_evp});
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
