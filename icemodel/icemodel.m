function [ice1, ice2, opts] = icemodel(opts)
   % ICEMODEL Simulate the phase change process in glacier ice.
   %
   % This function models the phase change process in melting glacier ice. It
   % uses iterative processes to update the temperature, liquid and ice fraction
   % at each time step. This model considers both the surface and subsurface
   % energy balance to simulate the phase change.
   %
   % Syntax:
   % [ice1, ice2] = ICEMODEL(opts)
   % [ice1, ice2, opts] = ICEMODEL(opts)
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
   % opts  - Finalized runtime configuration after icemodel.configureRun() has
   %         applied the last non-negotiable pre-execution updates.
   %
   % See also: skinmodel, icemodel.setopts
   %
   %#codegen

   %% INITIALIZE THE MODEL

   % Runtime configuration
   assertF on
   opts = icemodel.configureRun(opts);
   opts = icemodel.prepareRunOutput(opts);

   % LOAD PHYSICAL CONSTANTS AND PARAMETERS
   [cv_air, cv_ice, cv_liq, SB, k_liq, Lf, Ls, Lv, ro_air, ...
      ro_ice, ro_liq, roLf, roLs, roLv, Rv, Tf] ...
      = icemodel.physicalConstant( ...
      'cv_air','cv_ice','cv_liq','SB','k_liq','Lf','Ls','Lv','ro_air', ...
      'ro_ice','ro_liq','roLf','roLs','roLv','Rv','Tf');
   [emiss, fcp] = icemodel.parameterLookup('emiss', 'fcp');
   TINY = 1e-8;

   % UNPACK SOLVER OPTS
   [solver, seb_solver, maxiter, tol, alpha, use_aitken, jumpmax, ...
      cpl_maxiter, cpl_Ts_tol, cpl_seb_tol, cpl_alpha, cpl_aitken, cpl_jumpmax] ...
      = icemodel.getopts(opts, ...
      'solver', 'seb_solver', 'maxiter', 'tol', 'alpha', 'use_aitken', ...
      'jumpmax', 'cpl_maxiter', 'cpl_Ts_tol', 'cpl_seb_tol', 'cpl_alpha', ...
      'cpl_aitken', 'cpl_jumpmax');

   % LOAD THE FORCING DATA
   [tair, swd, lwd, albedo, wspd, rh, psfc, ppt, tppt, De, scoef, time] ...
      = METINIT(opts);

   % INITIALIZE THE SPECTRAL MODEL
   [I0, dz_spect, z_nodes_spect, z_edges_spect, tau_N, tau_S, solar_dwavel, ...
      k_bulk_lookup, r_eff] = EXTCOEFSINIT(opts, ro_ice);

   % INITIALIZE THE THERMAL MODEL
   [ice1, ice2, T_sfc, T_ice, f_ice, f_liq, r_eff, k_eff, fn, dz, delz, ...
      z_nodes, roL, liqflag, JJ, ~, Sp, Fc, Fp, TL, TH, f_ell_min, f_ell_max, ...
      f_ice_min, f_liq_res] = ICEINIT(opts, tair, r_eff);

   % INITIALIZE TIMESTEPPING
   [metstep, substep, numsteps, maxsubstep, dt, dt_FULL_STEP, ...
      numyears, numspinup] = INITTIMESTEPS(opts, time);
   if ~opts.saveflag && (numyears - numspinup) > 1
      ice1_all = [];
      ice2_all = [];
   end

   % INITIALIZE PAST VALUES
   [xT_sfc, xT_ice, xf_ice, xf_liq] = RESETSUBSTEP(T_sfc, T_ice, f_ice, f_liq);

   %% START TIMESTEPS OVER YEARS
   for thisyear = 1:numyears

      for timestep = 1:numsteps

         % INITIALIZE NEW TIMESTEP
         [dt_sum, n_subfail, ok_seb, ok_ieb, d_liq, d_evp, d_lyr] ...
            = NEWTIMESTEP(f_liq, solver);

         % SUBSURFACE SOLAR RADIATION SOURCE-TERM
         [Sc, chi] = SPECTRALSOURCETERM(swd(metstep), ...
            albedo(metstep), I0, dz_spect, tau_N, tau_S, solar_dwavel, ...
            dz, ro_ice * f_ice + ro_liq * f_liq, z_nodes, z_nodes_spect, ...
            z_edges_spect, k_bulk_lookup);

         % SURFACE TERMS (atmospheric vapor pressure fixed over this full step)
         ea = icemodel.surface.atmospheric_vapor_pressure(tair(metstep), ...
            rh(metstep), liqflag);
         snow_depth = 0.0;

         while dt_sum + TINY < dt_FULL_STEP
            
            % Update surface density for the surface turbulent heat flux scheme.
            ro_sfc = icemodel.surface.surface_bulk_density(f_ice(1), f_liq(1));

            if solver == 1 % Dirichlet iterated lagged closure mode

               % SURFACE TEMPERATURE
               k_eff = BULKTHERMALK(T_ice, f_ice, f_liq, ro_ice, k_liq);
               [T_sfc, ok_seb] = SEBSOLVE(T_sfc, tair(metstep), swd(metstep), ...
                  lwd(metstep), albedo(metstep), wspd(metstep), ppt(metstep), ...
                  tppt(metstep), psfc(metstep), De(metstep), ea, chi, roL, ...
                  scoef, liqflag, T_ice, k_eff, dz, ro_sfc, snow_depth, ...
                  seb_solver, opts.debug, opts);

               % SUBSURFACE ENERGY BALANCE
               [T_ice, f_ice, f_liq, k_eff, ok_ieb, n_iters] ...
                  = ICEENBAL(T_sfc, T_ice, f_ice, f_liq, Fc, Fp, Sc, Sp, ...
                  dz, delz, fn, dt, JJ, k_liq, cv_ice, cv_liq, ro_ice, ...
                  ro_liq, Ls, Lf, roLf, Tf, fcp, TL, TH, f_ell_min, ...
                  f_ell_max, solver, tol, maxiter, use_aitken, jumpmax, ...
                  opts.debug);

            elseif solver > 1 % Robin single-sweep and strong coupling modes

               % COUPLED SURFACE-SUBSURFACE ENERGY BALANCE
               [T_sfc, T_ice, f_ice, f_liq, k_eff, ok_ieb, n_iters] ...
                  = ICEEBSOLVE(T_ice, f_ice, f_liq, dz, delz, fn, Sc, dt, JJ, T_sfc, ...
                  k_liq, cv_ice, cv_liq, ro_ice, ro_liq, Ls, Lf, roLf, Rv, ...
                  Tf, fcp, TL, TH, f_ell_min, f_ell_max, tair(metstep), ...
                  swd(metstep), lwd(metstep), albedo(metstep), wspd(metstep), ...
                  ppt(metstep), tppt(metstep), psfc(metstep), De(metstep), ...
                  ea, cv_air, emiss, SB, roL, scoef, chi, liqflag, solver, ...
                  tol, maxiter, alpha, use_aitken, jumpmax, cpl_Ts_tol, ...
                  cpl_seb_tol, cpl_maxiter, cpl_alpha, cpl_aitken, cpl_jumpmax, ...
                  opts.debug, ro_sfc, snow_depth, opts);
            end
            ok = ok_seb && ok_ieb;

            % CHECK SUBSTEP FAILURE (shorten dt and restart substep on failure)
            [T_sfc, T_ice, f_ice, f_liq, n_subfail, substep, dt, ok] ...
               = CHECKSUBSTEP(T_sfc, T_ice, f_ice, f_liq, xT_sfc, xT_ice, xf_ice, xf_liq, ...
               ro_ice, ro_liq, dt_sum, dt, dt_FULL_STEP, timestep, numsteps, ...
               substep, maxsubstep, n_subfail, opts.debug, eps, ok);

            if ~ok
               continue
            end

            % UPDATE POTENTIAL SURFACE NET VAPOR FLUX
            ro_sfc = icemodel.surface.surface_bulk_density(f_ice(1), f_liq(1));
            [d_pevp, ~, ~, ~] = icemodel.surface.potential_surface_vapor_tendency( ...
               T_sfc, Tf, tair(metstep), wspd(metstep), psfc(metstep), ...
               De(metstep), ea, cv_air, roL, scoef, liqflag, ro_sfc, ...
               snow_depth, opts, Lv, ro_liq, dt, dz(1));

            % UPDATE MASS BALANCE FLUXES
            [T_ice, f_ice, f_liq, d_liq, d_evp, d_lyr] = ICEMF(T_ice, f_ice, f_liq, ...
               ro_ice, ro_liq, cv_ice, cv_liq, Lf, Ls, Lv, Tf, TL, fcp, ...
               xf_liq, Sc, Sp, JJ, f_ice_min, dz(1), d_pevp, d_liq, d_evp, ...
               d_lyr, f_ell_min, f_liq_res);

            % CHECKPOINT STATE AND SUBSTEP TIME
            [xT_sfc, xT_ice, xf_ice, xf_liq, dt_sum, dt, liqflag, roL] ...
               = UPDATESUBSTEP(T_sfc, T_ice, f_ice, f_liq, dt_FULL_STEP, dt_sum, ...
               dt, TINY, ro_ice, ro_liq, ro_air, cv_ice, cv_liq, roLv, roLs);
         end

         % Error if dt accumulation exceeds full step
         assertF(@() dt_sum < dt_FULL_STEP + 2 * TINY)

         % UPDATE GRAIN SIZE VIA VAPOR MASS TRANSFER
         r_eff = VAPORTRANSFER(T_ice, T_sfc, f_ice, f_liq, r_eff, dz, delz, ...
            fn, dt_FULL_STEP);

         % DIAGNOSE SURFACE FLUXES
         ro_sfc = icemodel.surface.surface_bulk_density(f_ice(1), f_liq(1));
         [Qe, Qh, Qc, Qm, ~, Qbal] = SEBFLUX(T_ice, T_sfc, tair(metstep), ...
            swd(metstep), lwd(metstep), albedo(metstep), wspd(metstep), ...
            ppt(metstep), tppt(metstep), psfc(metstep), De(metstep), ea, ...
            Tf, k_eff, dz, roL, chi, scoef, liqflag, ro_sfc, snow_depth, ...
            opts);

         % SAVE OUTPUT IF SPINUP IS FINISHED
         if thisyear > numspinup

            if strcmp(opts.output_profile, 'minimal')

               [ice1, ice2] = SAVEOUTPUT(timestep, ice1, ice2, ...
                  opts.vars1, opts.vars2, ...
                  {T_sfc},  ...
                  {T_ice, f_ice, f_liq, d_liq, d_evp});

            elseif strcmp(opts.output_profile, 'standard')

               [ice1, ice2] = SAVEOUTPUT(timestep, ice1, ice2, ...
                  opts.vars1, opts.vars2, ...
                  {T_sfc, Qm, Qe, Qh, Qc, chi, Qbal, dt_sum, ok_seb, ok_ieb, n_iters}, ...
                  {T_ice, f_ice, f_liq, d_liq, d_evp, d_lyr, Sc, r_eff});
            end
         end

         % MOVE TO THE NEXT TIMESTEP
         [metstep, substep, dt] = NEXTSTEP(metstep, substep, dt_FULL_STEP, ...
            maxsubstep, ok, n_subfail, n_iters);

      end % timesteps (one year)

      if isfield(opts, 'saverestart') && opts.saverestart
         icemodel.saveRestartState(opts, opts.simyears(thisyear), ...
            T_ice, f_ice, f_liq, T_sfc, r_eff);
      end

      % RESTART THE MET DATA STEP INDEX DURING SPIN UP
      if thisyear <= numspinup
         continue
      end

      % Concatenate yearly raw output when running multi-year simulations
      % without writing each year to disk.
      if ~opts.saveflag && numyears - numspinup > 1
         [ice1_all, ice2_all] = icemodel.concatoutput(ice1_all, ice2_all, ...
            ice1, ice2);
      end

      % WRITE TO DISK
      yridx = (thisyear-1)*numsteps+1:thisyear*numsteps;
      WRITEOUTPUT(ice1, ice2, opts, thisyear, ...
         time(yridx), swd(yridx), lwd(yridx), albedo(yridx))
   end

   if ~opts.saveflag && numyears - numspinup > 1
      ice1 = ice1_all;
      ice2 = ice2_all;
   end
end
