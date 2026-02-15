function [ice1, ice2] = icemodel(opts)
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
   %        * smbmodel         - Simulation model identifier
   %        * simyears         - Years for which the simulation is done
   %        * forcings         - Type of forcing data used
   %        * userdata         - User-defined data type
   %        * uservars         - User-defined variables
   %        * saveflag         - Flag indicating if data should be saved
   %        * testname         - Name of the test (default: 'none')
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
   %
   %%#codegen

   %% INITIALIZE THE MODEL

   debug = true;
   assertF on

   % LOAD PHYSICAL CONSTANTS AND PARAMETERS
   [cv_air, cv_ice, cv_liq, emiss, SB, epsilon, fcp, k_liq, Lf, ...
      Ls, Lv, ro_air, ro_ice, ro_liq, roLf, roLs, roLv, Rv, Tf] ...
      = icemodel.physicalConstant( ...
      'cv_air','cv_ice','cv_liq','emiss', 'SB','epsilon','fcp','k_liq','Lf', ...
      'Ls','Lv','ro_air','ro_ice','ro_liq','roLf','roLs','roLv', 'Rv','Tf');
   TINY = 1e-8;

   % LOAD THE FORCING DATA
   [tair, swd, lwd, albedo, wspd, rh, psfc, ppt, tppt, De, scoef, time] ...
      = METINIT(opts, 1);

   % INITIALIZE THE THERMAL MODEL
   [ice1, ice2, T, f_ice, f_liq, k_eff, fn, dz, delz, roL, liqflag, Ts, JJ, ...
      Sc, Sp, Fc, Fp, TL, TH, f_ell_min, f_ell_max, f_ice_min, f_liq_res] ...
      = ICEINIT(opts, tair);

   % INITIALIZE THE SPECTRAL MODEL
   [Q0, dz_spect, spect_N, spect_S, solardwavl] = EXTCOEFSINIT(opts, ro_ice);

   % INITIALIZE TIMESTEPPING
   [metiter, subiter, maxiter, maxsubiter, dt, dt_FULL_STEP, ...
      numyears, numspinup] = INITTIMESTEPS(opts, time);

   % INITIALIZE PAST VALUES
   [xTs, xT, xf_ice, xf_liq] = RESETSUBSTEP(Ts, T, f_ice, f_liq);

   bc = opts.bc_type;
   ok = true;

   %% START ITERATIONS OVER YEARS
   for thisyear = 1:numyears

      for iter = 1:maxiter

         % INITIALIZE NEW TIMESTEP
         [dt_sum, subfail, OK, d_liq, d_evp, d_lyr] = NEWTIMESTEP(f_liq);

         % SUBSURFACE SOLAR RADIATION SOURCE-TERM
         [Sc, chi] = UPDATEEXTCOEFS(swd(metiter), albedo(metiter), ...
            Q0, dz_spect, spect_N, spect_S, solardwavl, Sc, dz, ...
            ro_ice * f_ice + ro_liq * f_liq);

         % SURFACE TEMPERATURE
         ea = VAPPRESS(tair(metiter), Tf, liqflag) * rh(metiter) / 100;
         if bc == 1
            [Ts, ok] = SEBSOLVE(tair(metiter), swd(metiter), lwd(metiter), ...
               albedo(metiter), wspd(metiter), ppt(metiter), tppt(metiter), ...
               psfc(metiter), De(metiter), ea, cv_air, cv_liq, emiss, SB, Tf, ...
               chi, roL, scoef, liqflag, Ts, T, k_eff, dz, opts.seb_solver);

         elseif bc == 2
            [Fc, Fp] = SFCFLIN(tair(metiter), swd(metiter), lwd(metiter), ...
               albedo(metiter), wspd(metiter), psfc(metiter), De(metiter), ...
               ea, cv_air, emiss, SB, roL, scoef, chi, Tf, Ts, liqflag);
         end

         while dt_sum + TINY < dt_FULL_STEP

            % SUBSURFACE ENERGY BALANCE
            [T, f_ice, f_liq, k_eff, OK, N, a1] = ICEENBAL(T, f_ice, f_liq, ...
               dz, delz, fn, Sc, dt, JJ, Ts, k_liq, cv_ice, cv_liq, ro_ice, ...
               ro_liq, Ls, Lf, roLf, Rv, Tf, fcp, TL, TH, f_ell_min, ...
               f_ell_max, Fc, Fp, bc);

            if debug == true
               % PROGRESS MESSAGE (SLOWS DOWN THE CODE A LOT)
               fprintf('iter = %d (%.2f%%), dt = %.0f, success = %s\n', ...
                  iter, 100*iter/maxiter, dt, mat2str(all(OK)))
               if OK
                  assertF(@() all(f_ice + f_liq * ro_liq / ro_ice <= 1 + eps))
               end
            end

            % PHASE BOUNDARY OVERSHOOT, DECREASE THE TIME STEP AND START OVER
            if not(OK)
               [Ts, T, f_ice, f_liq, subfail, subiter, dt] ...
                  = RESETSUBSTEP(xTs, xT, xf_ice, xf_liq, dt_FULL_STEP, ...
                  subiter, maxsubiter, subfail, dt_sum);
               if subfail < maxsubiter
                  continue
               end
            end

            % UPDATE SEB LINEARIZATION (substep update for lagged-robin bc)
            if bc == 2
               Ts = (Fc + a1 * T(1)) / (a1 - Fp);
               [Fc, Fp] = SFCFLIN(tair(metiter), swd(metiter), lwd(metiter), ...
                  albedo(metiter), wspd(metiter), psfc(metiter), De(metiter), ...
                  ea, cv_air, emiss, SB, roL, scoef, chi, Tf, Ts, liqflag);
            end

            % UPDATE SURFACE FLUXES
            [Qe, Qh, Qc, Qm, ~, balance] = SEBFLUX(T, Ts, tair(metiter), ...
               swd(metiter), lwd(metiter), albedo(metiter), wspd(metiter), ...
               ppt(metiter), tppt(metiter), psfc(metiter), De(metiter), ea, ...
               Tf, k_eff, dz, cv_air, cv_liq, roL, emiss, SB, chi, epsilon, ...
               scoef, liqflag);

            % COMPUTE MASS BALANCE
            [T, f_ice, f_liq, d_liq, d_evp, d_lyr] = ICEMF(T, f_ice, f_liq, ...
               ro_ice, ro_liq, cv_ice, cv_liq, Lf, Ls, Lv, Tf, TL, fcp, ...
               xf_liq, Sc, Sp, JJ, f_ice_min, dz(1), dt, Qe, d_liq, d_evp, ...
               d_lyr, f_ell_min, f_liq_res);

            % UPDATE DENSITY, HEAT CAPACITY, AND SUBSTEP TIME
            [xTs, xT, xf_ice, xf_liq, dt_sum, dt, liqflag, roL] ...
               = UPDATESUBSTEP(Ts, T, f_ice, f_liq, dt_FULL_STEP, dt_sum, ...
               dt, TINY, ro_ice, ro_liq, ro_air, cv_ice, cv_liq, roLv, roLs);
         end

         assertF(@() dt_sum < dt_FULL_STEP + 2 * TINY)

         % SAVE OUTPUT IF SPINUP IS FINISHED
         if thisyear >= numspinup

            if strcmp('sector', opts.sitename)

               [ice1, ice2] = SAVEOUTPUT(iter, ice1, ice2, ...
                  opts.vars1, opts.vars2, ...
                  {Ts},  ...
                  {T, f_ice, f_liq, d_liq, d_evp});
            else

               [ice1, ice2] = SAVEOUTPUT(iter, ice1, ice2, ...
                  opts.vars1, opts.vars2, ...
                  {Ts, Qm, Qe, Qh, Qc, chi, balance, dt_sum, ok, OK, N}, ...
                  {T, f_ice, f_liq, d_liq, d_evp, d_lyr, Sc});
            end
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
