function [ice1, ice2, numfail] = icemodel(opts) %#codegen
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
   [tair, swd, lwd, albedo, wspd, rh, psfc, ppt, tppt, De, scoef, time] ...
      = METINIT(opts, 1);

   % INITIALIZE THE THERMAL MODEL
   [ice1, ice2, T, f_ice, f_liq, k_eff, fn, dz, delz, roL, liqflag, Ts, JJ, ...
      Sc, Sp, Fc, Fp, TL, TH, flmin, flmax, f_min, liqresid, ro_iwe, ro_wie] ...
      = ICEINIT(opts, tair);

   % INITIALIZE THE SPECTRAL MODEL
   [Q0, dz_spect, spect_N, spect_S, solardwavl] = EXTCOEFSINIT(opts, ro_ice);

   % INITIALIZE TIMESTEPPING
   [metiter, subiter, maxiter, maxsubiter, dt, dt_FULL_STEP, ...
      numyears, numspinup] = INITTIMESTEPS(opts, time);

   % INITIALIZE PAST VALUES
   [xT, xf_ice, xf_liq] = RESETSUBSTEP(T, f_ice, f_liq);

   bc = opts.bc_type;
   ok = true;
   % zD = dz(1);

   numfail = zeros(maxiter, 1);

   %% START ITERATIONS OVER YEARS
   for thisyear = 1:numyears

      for iter = 1:maxiter

         % INITIALIZE NEW TIMESTEP
         [dt_sum, subfail, OK, d_liq, d_drn, d_evp] = NEWTIMESTEP(f_liq);

         % SUBSURFACE SOLAR RADIATION SOURCE-TERM
         [Sc, chi] = UPDATEEXTCOEFS(swd(metiter), albedo(metiter), ...
            Q0, dz_spect, spect_N, spect_S, solardwavl, Sc, dz, ...
            ro_liq * f_liq + ro_ice * f_ice);

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
               ea, cv_air, emiss, SB, roL, scoef, chi, Tf, T(1), liqflag);
         end

         % Use this if computing the top layer interaction length scale
         % xT1 = T(1);

         while dt_sum + TINY < dt_FULL_STEP

            % subfail goes above maxsubiter when T(1) is on the TL limit I think
            % and cannot ge

            % SUBSURFACE ENERGY BALANCE
            [T, f_ice, f_liq, k_eff, OK, N, a1] = ICEENBAL(T, f_ice, f_liq, ...
               dz, delz, fn, Sc, dt, JJ, Ts, k_liq, cv_ice, cv_liq, ro_ice, ...
               ro_liq, Ls, Lf, roLf, Rv, Tf, fcp, TL, TH, flmin, flmax, ...
               ro_iwe, ro_wie, Fc, Fp, bc);

            if debug == true
               % PROGRESS MESSAGE (SLOWS DOWN THE CODE A LOT)
               fprintf('iter = %d (%.2f%%), dt = %.0f, success = %s\n', ...
                  iter, 100*iter/maxiter, dt, mat2str(all(OK)))
               assertF(@() all(f_ice + f_liq * ro_wie <= 1 + eps))
            end

            % PHASE BOUNDARY OVERSHOOT, DECREASE THE TIME STEP AND START OVER
            if not(OK)
               [T, f_ice, f_liq, subfail, subiter, dt] ...
                  = RESETSUBSTEP(xT, xf_ice, xf_liq, dt_FULL_STEP, subiter, ...
                  maxsubiter, subfail, dt_sum);
               if subfail < maxsubiter
                  continue
               end
            end

            % UPDATE SURFACE FLUXES
            if bc == 2
               Ts = (Fc + a1 * T(1)) / (a1 - Fp);
            end
            [Qe, Qh, Qc, Qm, ~, balance] = SEBFLUX(T, Ts, tair(metiter), ...
               swd(metiter), lwd(metiter), albedo(metiter), wspd(metiter), ...
               ppt(metiter), tppt(metiter), psfc(metiter), De(metiter), ea, ...
               Tf, k_eff, dz, cv_air, cv_liq, roL, emiss, SB, chi, epsilon, ...
               scoef, liqflag);

            % COMPUTE MASS BALANCE
            [T, f_ice, f_liq, d_liq, d_evp, d_drn] = ICEMF(T, f_ice, f_liq, ...
               ro_ice, ro_liq, cv_ice, cv_liq, Lf, Ls, Lv, Tf, TL, fcp, ...
               xf_liq, Sc, Sp, JJ, f_min, dz(1), dt, Qe, ro_iwe, ...
               d_liq, d_drn, d_evp, flmin, liqresid);

            % UPDATE DENSITY, HEAT CAPACITY, AND SUBSTEP TIME
            [xT, xf_ice, xf_liq, dt_sum, dt, liqflag, roL] ...
               = UPDATESUBSTEP(T, f_ice, f_liq, dt_FULL_STEP, dt_sum, ...
               dt, TINY, ro_ice, ro_liq, ro_air, cv_ice, cv_liq, roLv, roLs);
         end

         % % Update the top layer interaction length
         % if swd(metiter) > 1 && (T(1) - xT1) > TINY && (T(1) - min(Ts, Tf)) > TINY
         %    zD = sqrt((k_eff(1) * dt * (T(1) - min(Ts, Tf))) ...
         %       / (cv_ice * f_ice(1) + cv_liq * f_liq(1) * (T(1) - xT1)));
         %    % zD = sqrt(k_eff(1) * dt / (cv_ice * f_ice(1) + cv_liq * f_liq(1)));
         % end

         assertF(@() dt_sum < dt_FULL_STEP + 2 * TINY)

         % SAVE OUTPUT IF SPINUP IS FINISHED
         if thisyear >= numspinup

            if strcmp('sector', opts.sitename)

               [ice1, ice2] = SAVEOUTPUT(iter, ice1, ice2, ...
                  opts.vars1, opts.vars2, ...
                  {Ts},  ...
                  {T, f_ice, f_liq, d_liq, d_drn, d_evp});
            else

               [ice1, ice2] = SAVEOUTPUT(iter, ice1, ice2, ...
                  opts.vars1, opts.vars2, ...
                  {Ts, Qm, Qe, Qh, Qc, chi, balance, dt_sum, ok, OK, N}, ...
                  {T, f_ice, f_liq, d_liq, d_drn, d_evp, Sc});
            end
         end

         % MOVE TO THE NEXT TIMESTEP
         [metiter, subiter, dt] = NEXTSTEP(metiter, subiter, ...
            dt, dt_FULL_STEP, maxsubiter, OK);

         numfail(iter) = subfail;

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
