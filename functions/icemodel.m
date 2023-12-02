function [ice1, ice2] = icemodel(opts) %#codegen
   % ICEMODEL Simulate the phase change process in melting glacier ice.
   %
   % This function models the phase change process in melting glacier ice. It
   % uses iterative processes to update the temperature, liquid and ice fraction
   % at each time step. This model considers both the surface and subsurface
   % energy balance to simulate the phase change.
   %
   % Syntax:
   % [ice1, ice2, met, opts] = ICEMODEL(opts)
   %
   % Inputs:
   % opts - A structure containing model options and parameters. Defined by the
   %        icemodel_opts function. The fields include:
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
   % opts  - Updated struct containing model options and parameters.
   %
   % See also: ICEMODEL_OPTS

   %% INITIALIZE THE MODEL

   assertF on
   % assertF(@() all(f_ice + f_liq * ro_wie <= 1))

   % Define a logical flag that sets the simulation model
   % isicemodel = strcmp(opts.simmodel, 'icemodel');

   % Load the physical constants and parameters
   [cv_air, cv_liq, cv_ice, emiss, SB, epsilon, fcp, k_liq, Lf, ...
      Ls, Lv, ro_air, ro_ice, ro_liq, roLf, roLs, roLv, Rv, Tf] ...
      = icemodel.physicalConstant( ...
      'cv_air','cv_liq','cv_ice','emiss', 'SB','epsilon','fcp','k_liq','Lf', ...
      'Ls','Lv','ro_air','ro_ice','ro_liq','roLf','roLs','roLv', 'Rv','Tf');
   TINY = 1e-8;

   % Load the forcing data
   [tair, swd, lwd, albedo, wspd, rh, psfc, De, time] = METINIT(opts, 1);

   % Initialize the ice column
   [f_ice, f_liq, T, TL, TH, flmin, flmax, ~, k_eff, dz, fn, delz, z_therm, ...
      dz_therm, dz_spect, JJ, ~, ~, Sp, scoef, ~, ro_iwe, ro_wie, xTsfc, ...
      xf_liq, roL, Qc, f_min, fopts, liqresid, liqflag, ice1, ice2] ...
      = ICEINIT(opts, tair);

   % INITIALIZE THE EXTINCTION COEFFICIENTS
   [I0, z_spect, spect_N, spect_S, solardwavl] = EXTCOEFSINIT(opts, ro_ice);

   % INITIALIZE TIMESTEPPING
   [metiter, subiter, maxiter, maxsubiter, dt_FULL_STEP, dt_min, dt_max, ...
      dt_new, numyears, ~, numspinup] = INITTIMESTEPS(opts, time);

   Tsfc = tair(1);
   flag = false(maxiter, 1);
   numiter = zeros(maxiter, 1);

   %% START ITERATIONS OVER YEARS
   for thisyear = 1:numyears

      iter = 1;

      while iter <= maxiter

         % INITIALIZE NEW SUBSTEP
         [dt_sum, subfail, OK, d_liq, d_drn, d_evp] = INITSUBSTEP(f_liq);

         % update the upper layer diffusion length scale
         % if sqrt(k_eff(1)*dt/(ro_sno(1)*cp_sno(1))) > dz(1)
         %    pause;
         % end

         % SUBSURFACE SOLAR RADIATION SOURCE-TERM
         if swd(metiter) > 0 % && isicemodel
            [Sc, chi] = UPDATEEXTCOEFS(swd(metiter), albedo(metiter), I0, ...
               dz, z_spect, dz_spect, z_therm, dz_therm, spect_N, spect_S, ...
               solardwavl, ro_liq * f_liq + ro_ice * f_ice);
         else
            Sc = 0.0 * dz;
            chi = 1.0;
         end

         % icemodel.verifyShortwaveBalance(Qsi, Sc, albedo, chi, dz, 1);

         % SURFACE TEMPERATURE
         ea = VAPPRESS(tair(metiter), Tf, liqflag) * rh(metiter) / 100;
         Ts = SEBSOLVE(tair(metiter), swd(metiter), lwd(metiter), ...
            albedo(metiter), wspd(metiter), psfc(metiter), De(metiter), ...
            ea, cv_air, emiss, SB, Tf, chi, roL, scoef, liqflag, ...
            Ts, T, k_eff, dz, opts.sebsolver);

            assertF(@() all(f_ice + f_liq * ro_wie <= 1))

            % SUBSURFACE ENERGY BALANCE
            [T, errH, errT, f_ice, f_liq, OK1] = ICEENBAL(T, f_ice, f_liq, ...
               k_liq, cv_ice, cv_liq, ro_ice, ro_liq, ro_air, Ls, Lf, roLf, ...
               Rv, Tf, dz, delz, fn, dt_new, JJ, Tsfc, Sc, fcp, TL, TH, ...
               flmin, flmax, ro_iwe, ro_wie);

            assertF(@() all(f_ice + f_liq * ro_wie <= 1))

            % ERROR MESSAGE (SLOWS DOWN THE CODE A LOT)
            % fprintf('iter = %d (%.2f%%), dt = %.0f, success = %s\n', ...
            %    iter,100*iter/maxiter,dt_new,mat2str(all(OK)))

            % PHASE BOUNDARY OVERSHOOT, DECREASE THE TIME STEP AND START OVER
            if not(OK) && subfail < maxsubiter

               [subfail, subiter, dt_new, T, Tsfc, f_ice, f_liq] ...
                  = RESETSUBSTEP( xT, xTsfc, xf_ice, xf_liq, dt_max, ...
                  subiter, maxsubiter, subfail);

               continue
            else

               % NOTE: if subiter == maxsubiter and OK is false, we have a
               % problem b/c T is likely complex and xT etc get set in
               % UPDATESUBSTEP, so either need to reduce the timestep further,
               % or add a check that only sets xT etc if OK is true, and
               % otherwise keeps state constant essentially skipping the step
               if subiter == maxsubiter || subiter == maxsubiter - 1
                  assertF(@() OK)
               end

               % UPDATE SURFACE FLUXES
               k_eff = GETGAMMA(T, f_ice, f_liq, ro_ice, k_liq, Ls, Rv, Tf);
               [Qe, Qh, Qc, Qm, Qf, balance] = SEBFLUX(T, Ts, tair(metiter), ...
                  swd(metiter), lwd(metiter), albedo(metiter), wspd(metiter), ...
                  psfc(metiter), De(metiter), ea, Tf, k_eff, dz, cv_air, roL, ...
                  emiss, SB, chi, epsilon, scoef, liqflag);

               % COMPUTE MASS BALANCE
               [T, f_ice, f_liq, d_liq, d_evp, d_drn, ~, lcflag] = ICEMF(T, f_ice, ...
                  f_liq, ro_ice, ro_liq, cv_ice, cv_liq, Lf, Ls, Lv, Tf, ...
                  TL, fcp, xf_liq, Sc, Sp, JJ, f_min, fopts, dz_therm, ...
                  dt_new, Qe, ro_iwe, d_liq, d_drn, d_evp, flmin, liqresid);

               assertF(@() all(f_ice + f_liq * ro_wie <= 1))

               % UPDATE DENSITY, HEAT CAPACITY, AND SUBSTEP TIME
               [xT, xTsfc, xf_ice, xf_liq, dt_sum, dt_new, liqflag, roL] ...
                  = UPDATESUBSTEP(T, Tsfc, f_ice, f_liq, dt_FULL_STEP, ...
                  dt_sum, dt_new, TINY, ro_ice, ro_liq, ro_air, cv_ice, ...
                  cv_liq, roLv, roLs);

               % zD = sqrt(k_eff(1)*dt/(ro_sno(1)*cp_sno(1)));
               % if zD > dz(1)
               % end
            end
         end

         assertF(@() dt_sum < dt_FULL_STEP + 2 * TINY)

         % SAVE OUTPUT IF SPINUP IS FINISHED
         if thisyear >= numspinup
            [ice1, ice2] = SAVEOUTPUT(iter, ice1, ice2, ...
               opts.vars1, opts.vars2, ...
               {Tsfc, Qm, Qf, Qe, Qh, Qc, chi, balance, dt_sum},  ...
               {T, f_ice, f_liq, d_liq, d_drn, d_evp, Sc, errH, errT, lcflag});

            % For a regional run:
            % {Tsfc},  ...
            %    {T, f_ice, f_liq, d_liq, d_drn, d_evp});
         end

         % MOVE TO THE NEXT TIMESTEP
         [iter, metiter, subiter, dt_new] = NEXTSTEP(iter, metiter, ...
            subiter, dt_new, dt_max, dt_min, OK);

         assertF(@() all(f_ice + f_liq * ro_wie <= 1))

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
