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

   debug = true;
   assertF on
   opts = icemodel.configureRun(opts);
   opts = icemodel.prepareRunOutput(opts);

   % LOAD PHYSICAL CONSTANTS AND PARAMETERS
   [cv_air, cv_ice, cv_liq, emiss, SB, epsilon, fcp, k_liq, Lf, ...
      Ls, Lv, ro_air, ro_ice, ro_liq, roLf, roLs, roLv, Rv, Tf] ...
      = icemodel.physicalConstant( ...
      'cv_air','cv_ice','cv_liq','emiss', 'SB','epsilon','fcp','k_liq','Lf', ...
      'Ls','Lv','ro_air','ro_ice','ro_liq','roLf','roLs','roLv', 'Rv','Tf');
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

   % INITIALIZE THE THERMAL MODEL
   [ice1, ice2, T, f_ice, f_liq, k_eff, fn, dz, delz, z_nodes, roL, ...
      liqflag, Ts, JJ, Sc, Sp, Fc, Fp, TL, TH, f_ell_min, f_ell_max, ...
      f_ice_min, f_liq_res] = ICEINIT(opts, tair);

   % INITIALIZE THE SPECTRAL MODEL
   [I0, dz_spect, z_nodes_spect, z_edges_spect, tau_N, tau_S, solar_dwavel] ...
      = EXTCOEFSINIT(opts, ro_ice);
   [spectral_variant, k_bulk_lookup] ...
      = configureSpectralVariant(opts, dz_spect, tau_N, tau_S, solar_dwavel);

   % INITIALIZE TIMESTEPPING
   [metstep, substep, numsteps, maxsubstep, dt, dt_FULL_STEP, ...
      numyears, numspinup] = INITTIMESTEPS(opts, time);
   if ~opts.saveflag && (numyears - numspinup) > 1
      ice1_all = [];
      ice2_all = [];
   end

   % INITIALIZE PAST VALUES
   [xTs, xT, xf_ice, xf_liq] = RESETSUBSTEP(Ts, T, f_ice, f_liq);

   %% START TIMESTEPS OVER YEARS
   for thisyear = 1:numyears

      for timestep = 1:numsteps

         % INITIALIZE NEW TIMESTEP
         [dt_sum, n_subfail, ok_seb, ok_ieb, d_liq, d_evp, d_lyr] ...
            = NEWTIMESTEP(f_liq, solver);

         % SUBSURFACE SOLAR RADIATION SOURCE-TERM
         switch spectral_variant
            case "inlined"
               [Sc, chi] = SPECTRALSOURCETERM_INLINE(swd(metstep), ...
                  albedo(metstep), I0, dz_spect, tau_N, tau_S, solar_dwavel, ...
                  dz, ro_ice * f_ice + ro_liq * f_liq, z_nodes, z_nodes_spect);
            case "functions"
               [Sc, chi] = SPECTRALSOURCETERM(swd(metstep), ...
                  albedo(metstep), I0, dz_spect, tau_N, tau_S, solar_dwavel, ...
                  dz, ro_ice * f_ice + ro_liq * f_liq, z_nodes, z_nodes_spect, ...
                  z_edges_spect, k_bulk_lookup);
            case "lookup"
               [Sc, chi] = SPECTRALSOURCETERM(swd(metstep), ...
                  albedo(metstep), I0, dz_spect, tau_N, tau_S, solar_dwavel, ...
                  dz, ro_ice * f_ice + ro_liq * f_liq, z_nodes, z_nodes_spect, ...
                  z_edges_spect, k_bulk_lookup);
            otherwise
               error('Unrecognized spectral variant')
         end

         % SURFACE TERMS (atmospheric vapor pressure fixed over this full step)
         ea = VAPPRESS(tair(metstep), Tf, liqflag) * rh(metstep) / 100;

         while dt_sum + TINY < dt_FULL_STEP

            if solver == 1 % Dirichlet iterated lagged closure mode

               % SURFACE TEMPERATURE
               k_eff = GETGAMMA(T, f_ice, f_liq, ro_ice, k_liq, Ls, Rv, Tf);
               [Ts, ok_seb] = SEBSOLVE(tair(metstep), swd(metstep), ...
                  lwd(metstep), albedo(metstep), wspd(metstep), ...
                  ppt(metstep), tppt(metstep), psfc(metstep), De(metstep), ...
                  ea, cv_air, cv_liq, emiss, SB, Tf, chi, roL, scoef, ...
                  liqflag, Ts, T, k_eff, dz, seb_solver);

               % SUBSURFACE ENERGY BALANCE
               [T, f_ice, f_liq, k_eff, ok_ieb, n_iters] ...
                  = ICEENBAL(T, f_ice, f_liq, dz, delz, fn, Sc, dt, JJ, Ts, ...
                  k_liq, cv_ice, cv_liq, ro_ice, ro_liq, Ls, Lf, roLf, Rv, ...
                  Tf, fcp, TL, TH, f_ell_min, f_ell_max, Fc, Fp, solver, ...
                  tol, maxiter, alpha, use_aitken, jumpmax);

            elseif solver > 1 % Robin single-sweep and strong coupling modes

               % COUPLED SURFACE-SUBSURFACE ENERGY BALANCE
               [Ts, T, f_ice, f_liq, k_eff, ok_ieb, n_iters] ...
                  = ICEEBSOLVE(T, f_ice, f_liq, dz, delz, fn, Sc, dt, JJ, Ts, ...
                  k_liq, cv_ice, cv_liq, ro_ice, ro_liq, Ls, Lf, roLf, Rv, ...
                  Tf, fcp, TL, TH, f_ell_min, f_ell_max, tair(metstep), ...
                  swd(metstep), lwd(metstep), albedo(metstep), wspd(metstep), ...
                  ppt(metstep), tppt(metstep), psfc(metstep), De(metstep), ...
                  ea, cv_air, emiss, SB, roL, scoef, chi, liqflag, solver, ...
                  tol, maxiter, alpha, use_aitken, jumpmax, cpl_Ts_tol, ...
                  cpl_seb_tol, cpl_maxiter, cpl_alpha, cpl_aitken, cpl_jumpmax);
            end
            ok = ok_seb && ok_ieb;

            % CHECK SUBSTEP FAILURE (shorten dt and restart substep on failure)
            [Ts, T, f_ice, f_liq, n_subfail, substep, dt, ok] ...
               = CHECKSUBSTEP(Ts, T, f_ice, f_liq, xTs, xT, xf_ice, xf_liq, ...
               ro_ice, ro_liq, dt_sum, dt, dt_FULL_STEP, timestep, numsteps, ...
               substep, maxsubstep, n_subfail, debug, eps, ok);

            if not(ok)
               continue
            end

            % UPDATE POTENTIAL SURFACE NET VAPOR FLUX
            Qe = LATENT(De(metstep), ...
               STABLEFN(tair(metstep), Ts, wspd(metstep), scoef), ...
               ea, VAPPRESS(Ts, Tf, liqflag), roL, epsilon, psfc(metstep));
            d_pevp = PEVAP(Qe, Lv, ro_liq, dt, dz(1));

            % UPDATE MASS BALANCE FLUXES
            [T, f_ice, f_liq, d_liq, d_evp, d_lyr] = ICEMF(T, f_ice, f_liq, ...
               ro_ice, ro_liq, cv_ice, cv_liq, Lf, Ls, Lv, Tf, TL, fcp, ...
               xf_liq, Sc, Sp, JJ, f_ice_min, dz(1), d_pevp, d_liq, d_evp, ...
               d_lyr, f_ell_min, f_liq_res);

            % CHECKPOINT STATE AND SUBSTEP TIME
            [xTs, xT, xf_ice, xf_liq, dt_sum, dt, liqflag, roL] ...
               = UPDATESUBSTEP(Ts, T, f_ice, f_liq, dt_FULL_STEP, dt_sum, ...
               dt, TINY, ro_ice, ro_liq, ro_air, cv_ice, cv_liq, roLv, roLs);
         end

         % Error if dt accumulation exceeds full step
         assertF(@() dt_sum < dt_FULL_STEP + 2 * TINY)

         % DIAGNOSE SURFACE FLUXES
         [Qe, Qh, Qc, Qm, ~, Qbal] = SEBFLUX(T, Ts, tair(metstep), ...
            swd(metstep), lwd(metstep), albedo(metstep), wspd(metstep), ...
            ppt(metstep), tppt(metstep), psfc(metstep), De(metstep), ea, ...
            Tf, k_eff, dz, cv_air, cv_liq, roL, emiss, SB, chi, epsilon, ...
            scoef, liqflag);

         % SAVE OUTPUT IF SPINUP IS FINISHED
         if thisyear > numspinup

            if strcmp(opts.output_profile, 'minimal')

               [ice1, ice2] = SAVEOUTPUT(timestep, ice1, ice2, ...
                  opts.vars1, opts.vars2, ...
                  {Ts},  ...
                  {T, f_ice, f_liq, d_liq, d_evp});

            elseif strcmp(opts.output_profile, 'standard')

               [ice1, ice2] = SAVEOUTPUT(timestep, ice1, ice2, ...
                  opts.vars1, opts.vars2, ...
                  {Ts, Qm, Qe, Qh, Qc, chi, Qbal, dt_sum, ok_seb, ok_ieb, n_iters}, ...
                  {T, f_ice, f_liq, d_liq, d_evp, d_lyr, Sc});
            end
         end

         % MOVE TO THE NEXT TIMESTEP
         [metstep, substep, dt] = NEXTSTEP(metstep, substep, ...
            dt, dt_FULL_STEP, maxsubstep, ok, n_subfail, n_iters);

      end % timesteps (one year)

      if isfield(opts, 'saverestart') && opts.saverestart
         icemodel.saveRestartState(opts, opts.simyears(thisyear), ...
            T, f_ice, f_liq, Ts);
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

function [variant, bulk_lookup] = configureSpectralVariant(opts, dz_spect, ...
      tau_N, tau_S, solar_dwavel)
   %CONFIGURESPECTRALVARIANT Resolve the spectral update implementation.

   % Use the lookup bulk-extinction path in production. Tests and perf reports
   % can set OPTS.test_spectral_variant to benchmark the exact organized path
   % or the preserved legacy inline path against that accepted default.
   if isfield(opts, 'test_spectral_variant')
      variant = lower(string(opts.test_spectral_variant));
   else
      variant = "lookup";
   end

   switch variant
      case {"lookup", ""}
         variant = "lookup";
         bulk_lookup = icemodel.makeBulkExtCoefsLookup(dz_spect, ...
            tau_N, tau_S, solar_dwavel);
      case "functions"
         bulk_lookup = struct([]);
      case "inlined"
         bulk_lookup = struct([]);
      otherwise
         error('unrecognized test_spectral_variant: %s', variant)
   end
end
