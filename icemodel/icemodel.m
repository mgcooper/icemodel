function [ice1, ice2, opts] = icemodel(opts)
   %ICEMODEL Simulate the phase change process in glacier ice.
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
   % opts  - Finalized runtime configuration after icemodel.configureRun()
   %         applies the required pre-execution updates.
   %
   % See also: skinmodel, icemodel.setopts
   %
   %#codegen

   %% INITIALIZE THE MODEL

   % Runtime configuration
   assertF on
   opts = icemodel.configureRun(opts);
   opts = icemodel.prepareRunOutput(opts);

   % Verification suite option to ask icemodel to return snow-model-like outputs
   % when no snow-model implementation is available.
   if isfield(opts, 'verification_synthetic_snow') ...
         && opts.verification_synthetic_snow
      [ice1, ice2, opts] = icemodel.verification.syntheticSnowModelRun(opts);
      return
   end

   % Option to log thf/seb diagnostics.
   use_thf_diag = strcmp(opts.output_profile, 'diagnostic');

   % INITIALIZE SOLVER SETTINGS
   % 'settings' holds the primary solver settings, 'settings0' holds the
   % defaults. checksubstep switches settings to recovery mode (underrelaxation
   % with acceleration off) after a coupler failure, and acceptsubstep restores
   % settings0 after a recovery.
   [settings, settings0] = icemodel.couplers.initialize_solver_settings(opts);

   % Unpack other model options
   f_ice_min = icemodel.getopts(opts, 'f_ice_min');
   TINY = 1e-8;

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
      f_ice, f_liq, ~, Sp, r_eff, k_eff, fn, dz, delz, z_nodes] ...
      = icemodel.column.initialize_column_state(opts, tair, r_eff);

   % INITIALIZE THE SURFACE STATE
   [ea_atm, ro_atm, cv_atm, nu_air, H_h, De_e, br_coefs] ...
      = icemodel.surface.initialize_surface_state(opts, tair, wspd, rh, psfc);

   % INITIALIZE TIMESTEPPING
   [metstep, substep, numsteps, dt, numyears, numspinup] ...
      = icemodel.timestepping.initialize_timesteps(opts, time);

   if ~opts.saveflag && (numyears - numspinup) > 1
      ice1_all = [];
      ice2_all = [];
   end

   % INITIALIZE PAST VALUES
   [xT_sfc, xT_ice, xf_ice, xf_liq, xk_eff] ...
      = icemodel.timestepping.resetsubstep(T_sfc, T_ice, f_ice, f_liq, k_eff);
   force_advance_streak_dt = 0.0;

   %% START TIMESTEPS OVER YEARS
   for thisyear = 1:numyears

      for timestep = 1:numsteps

         % INITIALIZE NEW TIMESTEP
         [dt_sum, d_liq, d_evp, d_lyr, d_rof, ...
            d_vap_liq, d_vap_ice, diag] ...
            = icemodel.timestepping.newtimestep(f_liq);

         % Zero the mass/energy budget for this forcing step.
         budget = icemodel.column.initialize_budget_state( ...
            T_ice, f_ice, f_liq, dz);

         % Get the forcing observation heights and corresponding
         % bulk-Richardson coefficients for this forcing step.
         step_opts = icemodel.surface.step_observation_heights(opts, metstep);
         br_coefs_step = br_coefs(min(metstep, size(br_coefs, 1)), :);

         % SUBSURFACE SOLAR RADIATION SOURCE-TERM
         [Sc, chi] = icemodel.column.shortwave_source_term( ...
            swd(metstep), albedo(metstep), I0, dz_spect, tau_N, tau_S, ...
            solar_dwavel, dz, z_nodes, z_nodes_spect, z_edges_spect, ...
            icemodel.column.bulk_density(f_ice, f_liq), k_bulk_lookup);

         snow_depth = icemodel.surface.resolve_forcing_snow_depth( ...
            forcing_snow_depth, metstep, opts.use_forcing_snow_depth_for_thf);

         while dt_sum + TINY < settings.dt_full_step

            % UPDATE SURFACE STATE FOR CURRENT FORCING AND COLUMN STATE
            [liqflag, ro_sfc, hv_atm, H_e, f_res_por] = ...
               icemodel.surface.update_surface_state( ...
               f_ice(1), f_liq(1), ro_atm(metstep), De_e(metstep), ...
               snow_depth, step_opts);

            if settings.solver <= 1

               % COUPLED DIRICHLET SURFACE-SUBSURFACE ENERGY BALANCE
               [T_sfc, T_ice, f_ice, f_liq, k_eff, U_vap, L_vap, ...
                  diag.substep] ...
                  = icemodel.couplers.solve_surface_column_dirichlet( ...
                  T_sfc, T_ice, f_ice, f_liq, Sc, Sp, dz, delz, fn, dt, ...
                  tair(metstep), swd(metstep), lwd(metstep), albedo(metstep), ...
                  wspd(metstep), ppt(metstep), tppt(metstep), psfc(metstep), ...
                  ea_atm(metstep), ro_atm(metstep), cv_atm(metstep), ...
                  nu_air(metstep), H_h(metstep), H_e, hv_atm, br_coefs_step, ...
                  liqflag, chi, ro_sfc, snow_depth, f_res_por, settings, ...
                  step_opts);

            elseif settings.solver > 1

               % COUPLED ROBIN SURFACE-SUBSURFACE ENERGY BALANCE
               [T_sfc, T_ice, f_ice, f_liq, k_eff, U_vap, L_vap, ...
                  diag.substep] ...
                  = icemodel.couplers.solve_surface_column_robin( ...
                  T_sfc, T_ice, f_ice, f_liq, Sc, Sp, dz, delz, fn, dt, ...
                  tair(metstep), swd(metstep), lwd(metstep), albedo(metstep), ...
                  wspd(metstep), ppt(metstep), tppt(metstep), psfc(metstep), ...
                  ea_atm(metstep), ro_atm(metstep), cv_atm(metstep), ...
                  nu_air(metstep), H_h(metstep), H_e, hv_atm, br_coefs_step, ...
                  liqflag, chi, ro_sfc, snow_depth, f_res_por, settings, ...
                  step_opts);
            end

            % CHECK SUBSTEP FAILURE
            % Retry the substep on failure using this decision tree:
            % On a coupler failure, retry once in recovery mode.
            % On any of the following, retry with shortened dt:
            % - a surface-solve failure;
            % - an inner-solve failure;
            % - a failed retry; or
            % - a coupler failure after a failed recovery-mode attempt.
            % After maxsubstep failures, force advance.
            [T_sfc, T_ice, f_ice, f_liq, k_eff, substep, dt, ...
               ok, forced_advance, force_advance_streak_dt, settings, ...
               diag] = icemodel.timestepping.checksubstep( ...
               T_sfc, T_ice, f_ice, f_liq, k_eff, xT_sfc, xT_ice, ...
               xf_ice, xf_liq, xk_eff, dt_sum, dt, timestep, ...
               numsteps, substep, force_advance_streak_dt, 'icemodel', ...
               settings, settings0, diag);

            if ~ok
               continue
            end

            % On forced advance, update time, checkpoint the accepted state that
            % checksubstep just restored, and skip downstream physics.
            if forced_advance
               [xT_sfc, xT_ice, xf_ice, xf_liq, xk_eff, dt_sum, dt, ...
                  settings, diag] = ...
                  icemodel.timestepping.acceptsubstep( ...
                  T_sfc, T_ice, f_ice, f_liq, k_eff, dt_sum, ...
                  dt, TINY, settings, settings0, diag);
               continue
            end

            % Checkpoint melt/freeze phase change before surface vapor exchange.
            budget = icemodel.column.accumulate_phase_budget(budget, ...
               xT_ice, xf_ice, xf_liq, T_ice, f_ice, f_liq, dz);

            % UPDATE POTENTIAL SURFACE NET VAPOR FLUX
            d_pevp = icemodel.surface.potential_surface_vapor_demand( ...
               T_sfc, tair(metstep), wspd(metstep), psfc(metstep), ...
               ea_atm(metstep), ro_atm(metstep), cv_atm(metstep), ...
               nu_air(metstep), H_h(metstep), H_e, hv_atm, br_coefs_step, ...
               liqflag, f_ice(1), f_liq(1), dt, dz(1), snow_depth, step_opts);

            % APPLY AND BUDGET THE SURFACE MASS BALANCE
            [T_ice, f_ice, f_liq, d_liq, d_evp, d_rof, d_vap_liq, ...
               d_vap_ice, d_vap, budget] ...
               = icemodel.column.budget_surface_mass_balance( ...
               T_ice, f_ice, f_liq, xf_liq, d_pevp, d_liq, d_evp, ...
               d_rof, d_vap_liq, d_vap_ice, f_res_por, f_ice_min, ...
               budget, dz);

            % MOVE VAPOR BETWEEN THE CELLS
            [f_ice, f_liq, d_vap_liq, d_vap_ice, budget] ...
               = icemodel.column.couple_vapor_step( ...
               f_ice, f_liq, U_vap, L_vap, d_vap_liq, d_vap_ice, ...
               dz, dt, f_ice_min, f_res_por, budget);

            % UPDATE GRAIN SIZE FROM THIS SUBSTEP'S VAPOR EXCHANGE
            r_eff = icemodel.column.update_grain_radius( ...
               r_eff, f_liq, U_vap, d_vap, dz(1), dt);

            % REMESH THIN LAYERS AFTER THE MASS-BALANCE UPDATE
            [T_ice, f_ice, f_liq, Sc, Sp, d_lyr, budget] ...
               = icemodel.column.merge_thin_layers(T_ice, f_ice, f_liq, ...
               Sc, Sp, dz(1), d_pevp, d_lyr, f_ice_min, budget);

            % Vapor exchange and remeshing change column state, so k_eff needs
            % to be updated for the checkpoint and the diagnosed conduction.
            k_eff = icemodel.column.bulk_thermal_conductivity( ...
               T_ice, f_ice, f_liq, 0);

            % CHECKPOINT STATE AND SUBSTEP TIME
            [xT_sfc, xT_ice, xf_ice, xf_liq, xk_eff, dt_sum, dt, ...
               settings, diag] = ...
               icemodel.timestepping.acceptsubstep(T_sfc, T_ice, f_ice, ...
               f_liq, k_eff, dt_sum, dt, TINY, settings, settings0, diag);
         end

         % Error if dt accumulation exceeds a full step.
         assertF(@() dt_sum < settings.dt_full_step + 2 * TINY)

         % Close the forcing-step budget after the substep loop.
         budget = icemodel.column.finalize_budget_state( ...
            budget, T_ice, f_ice, f_liq, dz);

         % DIAGNOSE SURFACE ENERGY BALANCE
         [Qe, Qh, Qc, Qsn, Qln, Qa, Qm, Qf, Qbal] ...
            = icemodel.surface.diagnose_surface_energy_balance(T_sfc, ...
            tair(metstep), swd(metstep), lwd(metstep), albedo(metstep), ...
            wspd(metstep), ppt(metstep), tppt(metstep), psfc(metstep), ...
            ea_atm(metstep), ro_atm(metstep), cv_atm(metstep), ...
            nu_air(metstep), H_h(metstep), H_e, hv_atm, br_coefs_step, ...
            liqflag, chi, T_ice, k_eff, dz, ro_sfc, snow_depth, step_opts);

         % Build the thf diagnostic struct if requested.
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

            % Build the same surface_state fields regardless of
            % opts.output_profile, so codegen compiles one fixed struct shape.
            surface_state = budget;
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
            surface_state.Tsfc_converged = diag.ok_seb;
            surface_state.Tice_converged = diag.ok_ieb;
            surface_state.Tice_numiter = diag.n_iters;
            surface_state.n_failed_substeps = diag.n_failed_substeps;
            surface_state.n_forced_advances = diag.n_forced_advances;
            surface_state.cpl_iters = diag.cpl_iters;
            surface_state.cpl_res = diag.cpl_res;
            surface_state.seb_res = diag.seb_res;
            surface_state.cpl_recovery_count = diag.cpl_recovery_count;
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
               'df_vap_liq', d_vap_liq, ...
               'df_vap_ice', d_vap_ice, ...
               'Sc', Sc, ...
               'r_eff', r_eff);

            % Assemble the outputs.
            [data1, data2] = icemodel.buildOutputPayload(opts, ...
               surface_state, subsurface_state, thf_diag);

            % Update the requested ice1 and ice2 outputs for this timestep.
            [ice1, ice2] = icemodel.updateoutput(timestep, ice1, ice2, ...
               opts.vars1, opts.vars2, data1, data2);
         end

         % MOVE TO THE NEXT TIMESTEP
         [metstep, substep, dt] = icemodel.timestepping.nexttimestep( ...
            metstep, substep, ok, settings, diag);

      end % timesteps (one year)

      if isfield(opts, 'saverestart') && opts.saverestart
         icemodel.saveRestartState(opts, opts.simyears(thisyear), ...
            T_ice, f_ice, f_liq, T_sfc, r_eff);
      end

      % RESTART THE MET DATA STEP INDEX DURING SPIN UP
      if thisyear <= numspinup
         continue
      end

      % Concatenate yearly output when a run spans multiple years.
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
