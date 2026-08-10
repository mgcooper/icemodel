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

   % Verification can ask icemodel to return snow-model-like outputs before
   % the production snow physics exists. The bypass is used only by the
   % verification namespace, so the normal solver path stays unchanged.
   if isfield(opts, 'verification_synthetic_snow') ...
         && opts.verification_synthetic_snow
      [ice1, ice2, opts] = icemodel.verification.syntheticSnowModelRun(opts);
      return
   end

   TINY = 1e-8;

   % Both extra ledgers are diagnostic-only, so build them only there. They
   % are named apart because they are separate diagnostics: narrowing one gate
   % must not empty the other's channels.
   use_diagnostic_profile = strcmp(opts.output_profile, 'diagnostic');
   use_mass_budget = use_diagnostic_profile;
   use_thf_diag = use_diagnostic_profile;

   % UNPACK SOLVER OPTS
   [solver, maxiter, tol, alpha, use_aitken, jumpmax, cpl_maxiter, ...
      cpl_Ts_tol, cpl_seb_tol, cpl_alpha, cpl_aitken, cpl_jumpmax, f_ice_min] ...
      = icemodel.getopts(opts, ...
      'solver', 'maxiter', 'tol', 'alpha', 'use_aitken', 'jumpmax', ...
      'cpl_maxiter', 'cpl_Ts_tol', 'cpl_seb_tol', 'cpl_alpha', 'cpl_aitken', ...
      'cpl_jumpmax', 'f_ice_min');

   % INITIALIZE THE FORCING DATA
   [tair, swd, lwd, albedo, wspd, ...
      rh, psfc, ppt, tppt, time, forcing_snow_depth, opts] ...
      = icemodel.surface.initialize_surface_forcings(opts);

   % INITIALIZE THE SPECTRAL MODEL
   [I0, dz_spect, z_nodes_spect, ...
      z_edges_spect, tau_N, tau_S, solar_dwavel, k_bulk_lookup, r_eff] ...
      = icemodel.radiation.initialize_spectral_model(opts);

   % INITIALIZE THE THERMAL MODEL
   [ice1, ice2, T_sfc, T_ice, ...
      f_ice, f_liq, ~, Sp, r_eff, k_eff, fn, dz, delz, z_nodes, ~] ...
      = icemodel.column.initialize_column_state(opts, tair, r_eff);

   % INITIALIZE THE SURFACE STATE
   [ea_atm, ro_atm, cv_atm, nu_air, H_h, De_e, br_coefs] ...
      = icemodel.surface.initialize_surface_state(opts, tair, wspd, rh, psfc);

   % INITIALIZE TIMESTEPPING
   [metstep, substep, numsteps, ...
      maxsubstep, dt, dt_FULL_STEP, numyears, numspinup] ...
      = icemodel.timestepping.initialize_timesteps(opts, time);

   if ~opts.saveflag && (numyears - numspinup) > 1
      ice1_all = [];
      ice2_all = [];
   end

   % INITIALIZE PAST VALUES
   [xT_sfc, xT_ice, xf_ice, xf_liq] ...
      = icemodel.timestepping.resetsubstep(T_sfc, T_ice, f_ice, f_liq);
   force_advance_streak_dt = 0.0;

   % Define the ledger and the handoffs the accumulators pass between each
   % other. They are written and read inside `if use_mass_budget` blocks, and
   % use_mass_budget is a runtime value, so they need a value on every path.
   [mass_energy_budget, solid_p, liquid_p, phase_solid, phase_liquid, ...
      vapor_solid, vapor_liquid] = icemodel.column.initialize_budget_state();

   %% START TIMESTEPS OVER YEARS
   for thisyear = 1:numyears

      for timestep = 1:numsteps

         % INITIALIZE NEW TIMESTEP
         [dt_sum, n_subfail, ok_seb, ok_ieb, d_liq, d_evp, d_lyr, d_rof] ...
            = icemodel.timestepping.newtimestep(f_liq, solver);

         % Zero the ledger for this forcing step and take its opening
         % storage. The channels accumulate across substeps, so they have to
         % start from zero again each step.
         if use_mass_budget
            mass_energy_budget = icemodel.column.initialize_budget_state();
            [mass_energy_budget.mass_budget_solid_start_mwe, ...
               mass_energy_budget.mass_budget_liquid_start_mwe] = ...
               icemodel.column.integrate_column_budget(T_ice, f_ice, f_liq, dz);
         end

         % Scalarize time-varying met observation heights and corresponding
         % bulk-Richardson coefficients prior to each forcing step.
         step_opts = icemodel.surface.step_observation_heights(opts, metstep);
         br_coefs_step = br_coefs(min(metstep, size(br_coefs, 1)), :);

         % SUBSURFACE SOLAR RADIATION SOURCE-TERM
         [Sc, chi] = icemodel.column.shortwave_source_term(swd(metstep), ...
            albedo(metstep), I0, dz_spect, tau_N, tau_S, solar_dwavel, ...
            dz, z_nodes, z_nodes_spect, z_edges_spect, ...
            icemodel.column.bulk_density(f_ice, f_liq), k_bulk_lookup);

         snow_depth = icemodel.surface.resolve_forcing_snow_depth( ...
            forcing_snow_depth, metstep, opts.use_forcing_snow_depth_for_thf);

         while dt_sum + TINY < dt_FULL_STEP

            % UPDATE SURFACE STATE FOR CURRENT FORCING AND COLUMN STATE
            [liqflag, ro_sfc, hv_atm, H_e, f_res_por] = ...
               icemodel.surface.update_surface_state( ...
               f_ice(1), f_liq(1), ro_atm(metstep), ...
               De_e(metstep), snow_depth, step_opts);

            if solver <= 1

               % COUPLED DIRICHLET SURFACE-SUBSURFACE ENERGY BALANCE
               [T_sfc, T_ice, f_ice, f_liq, k_eff, ok_seb, ok_ieb, ok_cpl, n_iters] ...
                  = icemodel.couplers.solve_surface_column_dirichlet( ...
                  T_sfc, T_ice, f_ice, f_liq, Sc, Sp, dz, delz, fn, dt, ...
                  tair(metstep), swd(metstep), lwd(metstep), albedo(metstep), ...
                  wspd(metstep), ppt(metstep), tppt(metstep), psfc(metstep), ...
                  ea_atm(metstep), ro_atm(metstep), cv_atm(metstep), ...
                  nu_air(metstep), H_h(metstep), H_e, hv_atm, br_coefs_step, ...
                  liqflag, chi, solver, tol, maxiter, alpha, use_aitken, ...
                  jumpmax, cpl_Ts_tol, cpl_seb_tol, cpl_maxiter, cpl_alpha, ...
                  cpl_aitken, cpl_jumpmax, ro_sfc, snow_depth, step_opts);

            elseif solver > 1

               % COUPLED ROBIN SURFACE-SUBSURFACE ENERGY BALANCE
               [T_sfc, T_ice, f_ice, f_liq, k_eff, ok_seb, ok_ieb, ok_cpl, n_iters] ...
                  = icemodel.couplers.solve_surface_column_robin( ...
                  T_sfc, T_ice, f_ice, f_liq, Sc, Sp, dz, delz, fn, dt, ...
                  tair(metstep), swd(metstep), lwd(metstep), albedo(metstep), ...
                  wspd(metstep), ppt(metstep), tppt(metstep), psfc(metstep), ...
                  ea_atm(metstep), ro_atm(metstep), cv_atm(metstep), ...
                  nu_air(metstep), H_h(metstep), H_e, hv_atm, br_coefs_step, ...
                  liqflag, chi, solver, tol, maxiter, alpha, use_aitken, ...
                  jumpmax, cpl_Ts_tol, cpl_seb_tol, cpl_maxiter, cpl_alpha, ...
                  cpl_aitken, cpl_jumpmax, ro_sfc, snow_depth, step_opts);
            end

            % Hitting max coupling iterations without ok_cpl is a substep fail.
            ok = ok_seb && ok_ieb && ok_cpl;

            % CHECK SUBSTEP FAILURE (shorten dt and restart substep on failure)
            [T_sfc, T_ice, f_ice, f_liq, n_subfail, substep, dt, ok, ...
               ~, force_advance_streak_dt] ...
               = icemodel.timestepping.checksubstep(T_sfc, T_ice, f_ice, ...
               f_liq, xT_sfc, xT_ice, xf_ice, xf_liq, dt_sum, dt, ...
               dt_FULL_STEP, timestep, numsteps, substep, maxsubstep, ...
               n_subfail, opts.debug, eps, ok, force_advance_streak_dt, ...
               dt_FULL_STEP, 'icemodel');

            if ~ok
               continue
            end

            % Checkpoint melt/freeze phase change before surface vapor exchange.
            if use_mass_budget
               [mass_energy_budget, solid_p, liquid_p, phase_solid, phase_liquid] ...
                  = icemodel.column.accumulate_phase_budget(mass_energy_budget, ...
                  xT_ice, xf_ice, xf_liq, T_ice, f_ice, f_liq, dz);
            end

            % UPDATE POTENTIAL SURFACE NET VAPOR FLUX
            [d_pevp, ~, ~, ~] ...
               = icemodel.surface.potential_surface_vapor_tendency( ...
               T_sfc, tair(metstep), wspd(metstep), psfc(metstep), ...
               ea_atm(metstep), ro_atm(metstep), cv_atm(metstep), ...
               nu_air(metstep), H_h(metstep), H_e, hv_atm, br_coefs_step, ...
               liqflag, f_ice(1), f_liq(1), dt, dz(1), snow_depth, step_opts);

            % UPDATE THE SURFACE MASS-BALANCE BUDGETS
            [T_ice, f_ice, f_liq, d_liq, d_evp, d_rof, d_sbl_err] ...
               = icemodel.column.budget_surface_mass_balance( ...
               T_ice, f_ice, f_liq, xf_liq, d_pevp, d_liq, d_evp, d_rof, ...
               f_res_por, f_ice_min);

            % Checkpoint realized vapor exchange and its input, overflow, and
            % signed unapplied energy (d_sbl_err).
            if use_mass_budget
               [mass_energy_budget, vapor_solid, vapor_liquid] ...
                  = icemodel.column.accumulate_vapor_budget(mass_energy_budget, ...
                  solid_p, liquid_p, T_ice, f_ice, f_liq, dz, ...
                  d_pevp, d_rof, d_sbl_err);
            end

            % REMESH THIN LAYERS AFTER THE MASS-BALANCE UPDATE
            % The eighth output builds the per-event remesh ledger. Requesting
            % it on every call costs about 6 us per substep, which measured as
            % +18% on a model year (63.7 s to 75.1 s, kanm 2016 solver 1), so
            % the seven-output form is used unless the ledger is wanted.
            if use_mass_budget
               [T_ice, f_ice, f_liq, Sc, Sp, d_lyr, ~, remesh] ...
                  = icemodel.column.merge_thin_layers( ...
                  T_ice, f_ice, f_liq, Sc, Sp, dz(1), d_pevp, d_lyr, ...
                  f_ice_min);

               % Numerical remeshing, domain exchange, and discrete grid
               % translation are tracked apart from the physical phase and
               % vapor increments.
               mass_energy_budget = icemodel.column.accumulate_remesh_budget( ...
                  mass_energy_budget, remesh, phase_solid + vapor_solid, ...
                  phase_liquid + vapor_liquid);
            else
               [T_ice, f_ice, f_liq, Sc, Sp, d_lyr, ~] ...
                  = icemodel.column.merge_thin_layers( ...
                  T_ice, f_ice, f_liq, Sc, Sp, dz(1), d_pevp, d_lyr, f_ice_min);
            end

            % CHECKPOINT STATE AND SUBSTEP TIME
            [xT_sfc, xT_ice, xf_ice, xf_liq, dt_sum, dt] = ...
               icemodel.timestepping.updatesubstep( ...
               T_sfc, T_ice, f_ice, f_liq, ...
               dt_FULL_STEP, dt_sum, dt, TINY);
         end

         % Error if dt accumulation exceeds full step
         assertF(@() dt_sum < dt_FULL_STEP + 2 * TINY)

         % Close the forcing-step storage endpoints after every accepted
         % substep, including any numerical remeshing/domain exchange.
         if use_mass_budget
            [mass_energy_budget.mass_budget_solid_end_mwe, ...
               mass_energy_budget.mass_budget_liquid_end_mwe] = ...
               icemodel.column.integrate_column_budget( ...
               T_ice, f_ice, f_liq, dz);
         end

         % UPDATE GRAIN SIZE VIA VAPOR MASS TRANSFER
         r_eff = icemodel.column.vapor_mass_transfer(T_ice, T_sfc, ...
            f_ice, f_liq, r_eff, dz, delz, fn, dt_FULL_STEP);

         % DIAGNOSE SURFACE ENERGY BALANCE
         [Qe, Qh, Qc, Qsn, Qln, Qa, Qm, Qf, Qbal] ...
            = icemodel.surface.diagnose_surface_energy_balance(T_sfc, ...
            tair(metstep), swd(metstep), lwd(metstep), albedo(metstep), ...
            wspd(metstep), ppt(metstep), tppt(metstep), psfc(metstep), ...
            ea_atm(metstep), ro_atm(metstep), cv_atm(metstep), ...
            nu_air(metstep), H_h(metstep), H_e, hv_atm, br_coefs_step, ...
            liqflag, chi, T_ice, k_eff, dz, ro_sfc, snow_depth, step_opts);

         % Build a detailed thf diagnostic ledger if requested.
         if use_thf_diag
            [~, ~, thf_diag] ...
               = icemodel.surface.diagnose_turbulent_heat_fluxes( ...
               icemodel.surface.physical_surface_temperature(T_sfc), ...
               tair(metstep), wspd(metstep), psfc(metstep), ...
               ea_atm(metstep), ro_atm(metstep), cv_atm(metstep), ...
               nu_air(metstep), H_h(metstep), H_e, hv_atm, br_coefs_step, ...
               liqflag, ro_sfc, snow_depth, step_opts);
         else
            thf_diag = struct([]);
         end

         % SAVE OUTPUT IF SPINUP IS FINISHED
         if thisyear > numspinup

            % Assemble one compile-time struct layout for every profile. The
            % output profile still selects the externally visible fields.
            % Branching here would leave surface_state without a single
            % inferable type for Coder, and would make a caller who names a
            % mass_budget_* channel in vars1 outside the diagnostic profile
            % hit unknownOutputField instead of reading the zeroed ledger.
            surface_state = mass_energy_budget;
            surface_state.Tsfc = T_sfc;
            surface_state.Qm = Qm;
            surface_state.Qf = Qf;
            surface_state.Qe = Qe;
            surface_state.Qh = Qh;
            surface_state.Qc = Qc;
            surface_state.Qsn = Qsn;
            surface_state.Qln = Qln;
            surface_state.Qa = Qa;
            surface_state.chi = chi;
            surface_state.balance = Qbal;
            surface_state.dt_sum = dt_sum;
            surface_state.Tsfc_converged = ok_seb;
            surface_state.Tice_converged = ok_ieb;
            surface_state.Tice_numiter = n_iters;
            surface_state.n_subfail = n_subfail;
            surface_state.ea_atm = ea_atm(metstep);
            surface_state.br_coefs_gamma = br_coefs_step(1);
            surface_state.br_coefs_b1_num = br_coefs_step(2);
            surface_state.br_coefs_b2_num = br_coefs_step(3);
            surface_state.hv_atm = hv_atm;
            surface_state.ro_sfc = ro_sfc;
            surface_state.df_rof = d_rof;

            subsurface_state = struct( ...
               'Tice', T_ice, ...
               'f_ice', f_ice, ...
               'f_liq', f_liq, ...
               'df_liq', d_liq, ...
               'df_evp', d_evp, ...
               'df_lyr', d_lyr, ...
               'Sc', Sc, ...
               'r_eff', r_eff);

            [data1, data2] = icemodel.buildOutputPayload(opts, ...
               surface_state, subsurface_state, thf_diag);

            [ice1, ice2] = icemodel.updateoutput(timestep, ice1, ice2, ...
               opts.vars1, opts.vars2, data1, data2);
         end

         % MOVE TO THE NEXT TIMESTEP
         [metstep, substep, dt] = icemodel.timestepping.nexttimestep( ...
            metstep, substep, dt_FULL_STEP, maxsubstep, ok, n_subfail, n_iters);

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

      icemodel.writeoutput(ice1, ice2, opts, thisyear, ...
         time(yridx), swd(yridx), lwd(yridx), albedo(yridx))
   end

   if ~opts.saveflag && numyears - numspinup > 1
      ice1 = ice1_all;
      ice2 = ice2_all;
   end
end
