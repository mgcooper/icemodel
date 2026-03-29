function report = plot_spectral_variant_profiles(kwargs)
   %PLOT_SPECTRAL_VARIANT_PROFILES Plot functions and lookup spectral profiles.
   %
   %  report = plot_spectral_variant_profiles()
   %  report = plot_spectral_variant_profiles(output_dir="/tmp/spectral_profiles")
   %
   % This tool builds one representative daytime spectral state, then plots:
   %  1. bulk extinction coefficients on the spectral grid
   %  2. net spectral flux on the spectral grid
   %  3. absorbed spectral divergence before thermal-grid collapse
   %  4. thermal-grid volumetric source term after collapse

   arguments

      kwargs.simyear (1, 1) double {mustBeInteger, mustBePositive} ...
         = 2016

      kwargs.output_dir (1, :) string ...
         = fullfile(icemodel.getpath('test'), ...
         'benchmarks', 'figures', 'spectral', 'manual')

      kwargs.smbmodel (1, :) string ...
         = "icemodel"

      kwargs.solver (1, 1) double {mustBeMember(kwargs.solver, [1 2 3])} ...
         = 3

      kwargs.metstep (1, 1) double {mustBeInteger, mustBePositive} ...
         = 49
   end

   % Install the suite config so the synthetic workspace uses the canonical
   % demo-backed test environment.
   [~, ~, ~, ~, suite_cleanup] = ...
      icemodel.test.helpers.bootstrapTestEnvironment(); %#ok<ASGLU>

   % Build one daytime synthetic state and the shared spectral geometry.
   workspace = icemodel.test.fixtures.makeSyntheticWorkspace(kwargs.simyear, ...
      configure=true, nsteps=96, dt_seconds=900);
   cleanup = onCleanup(@() icemodel.test.fixtures.cleanupSyntheticWorkspace( ...
      workspace));
   state = icemodel.test.fixtures.makeSyntheticColumnState(workspace, ...
      kwargs.smbmodel, solver=kwargs.solver, metstep=kwargs.metstep, ...
      include_spectral=true, testname='spectral_profile_plots');

   profiles = buildProfiles(state);

   % Ensure the output directory exists before exporting the figures.
   if exist(char(kwargs.output_dir), 'dir') ~= 7
      mkdir(char(kwargs.output_dir));
   end

   % Export linear-depth and log-depth views so the near-surface detail is
   % visible without giving up the full-column context.
   figs = struct();
   figs.bulkcoefs = fullfile(kwargs.output_dir, 'bulkcoefs.png');
   figs.bulkcoefs_logy = fullfile(kwargs.output_dir, 'bulkcoefs_logy.png');
   figs.netflux = fullfile(kwargs.output_dir, 'netflux_spectral.png');
   figs.netflux_logy = fullfile(kwargs.output_dir, 'netflux_spectral_logy.png');
   figs.divergence = fullfile(kwargs.output_dir, 'divergence_spectral.png');
   figs.divergence_logy = fullfile(kwargs.output_dir, ...
      'divergence_spectral_logy.png');
   figs.source = fullfile(kwargs.output_dir, 'source_therm.png');
   figs.source_logy = fullfile(kwargs.output_dir, 'source_therm_logy.png');

   exportTwoPanelPlot(profiles.z_nodes_spect, profiles.bulk_functions, ...
      profiles.bulk_lookup, figs.bulkcoefs, ...
      'Bulk Extinction Coefficients', '\eta', 'Spectral Cell Center z (m)');
   exportTwoPanelPlot(profiles.z_nodes_spect, profiles.bulk_functions, ...
      profiles.bulk_lookup, figs.bulkcoefs_logy, ...
      'Bulk Extinction Coefficients', '\eta', 'Spectral Cell Center z (m)', ...
      use_log_depth=true);
   exportTwoPanelPlot(profiles.z_flux_nodes, profiles.Q_functions, ...
      profiles.Q_lookup, figs.netflux, 'Net Spectral Flux', 'Q', ...
      'Spectral Interface z (m)');
   exportTwoPanelPlot(profiles.z_flux_nodes, profiles.Q_functions, ...
      profiles.Q_lookup, figs.netflux_logy, 'Net Spectral Flux', 'Q', ...
      'Spectral Interface z (m)', use_log_depth=true);
   exportTwoPanelPlot(profiles.z_nodes_spect, profiles.dQ_functions, ...
      profiles.dQ_lookup, figs.divergence, ...
      'Absorbed Spectral Divergence', '\DeltaQ', ...
      'Spectral Cell Center z (m)');
   exportTwoPanelPlot(profiles.z_nodes_spect, profiles.dQ_functions, ...
      profiles.dQ_lookup, figs.divergence_logy, ...
      'Absorbed Spectral Divergence', '\DeltaQ', ...
      'Spectral Cell Center z (m)', use_log_depth=true);
   exportTwoPanelPlot(profiles.z_nodes_therm, profiles.Sc_functions, ...
      profiles.Sc_lookup, figs.source, ...
      'Thermal-Grid Source Term', 'S_c', 'Therm Cell Center z (m)');
   exportTwoPanelPlot(profiles.z_nodes_therm, profiles.Sc_functions, ...
      profiles.Sc_lookup, figs.source_logy, ...
      'Thermal-Grid Source Term', 'S_c', 'Therm Cell Center z (m)', ...
      use_log_depth=true);

   report = struct();
   report.profiles = profiles;
   report.figures = figs;
   report.note = [ ...
      "The functions path is the exact organized production path. ", ...
      "The lookup path is plotted against it because the inlined legacy ", ...
      "path is visually indistinguishable from the exact functions path."];
   clear cleanup
end

function profiles = buildProfiles(s)
   %BUILDPROFILES Compute exact and lookup spectral profiles on one state.

   % Build the density profiles seen by the exact and lookup paths.
   ro_sno = s.ro_ice * s.f_ice + s.ro_liq * s.f_liq + ...
      s.ro_air * (1 - s.f_ice - s.f_liq);
   ro_sno_spect = max(interp1(s.z_nodes, ro_sno, s.z_nodes_spect, ...
      'nearest', 'extrap'), 300);

   % Build the k_bulk lookup table and an empty one for the exact path.
   k_lookup = icemodel.makeBulkExtCoefsLookup( ...
      s.dz_spect, s.tau_N, s.tau_S, s.solar_dwavel);
   k_lookup_empty = struct([]);

   % Build the functions-path and lookup bulk-extinction profiles.
   bulk_functions = BULKEXTCOEFS( ...
      s.dz_spect, ro_sno_spect, s.tau_N, s.tau_S, s.solar_dwavel);
   bulk_lookup = BULKEXTCOEFSLOOKUP(ro_sno_spect, k_lookup);

   % Solve the two-stream system for the exact and lookup bulk coefficients.
   Q_functions = SOLVETWOSTREAM( ...
      s.I0, s.albedo, bulk_functions, s.z_edges_spect);
   Q_lookup = SOLVETWOSTREAM( ...
      s.I0, s.albedo, bulk_lookup, s.z_edges_spect);
   dQ_functions = Q_functions(1:end - 1) - Q_functions(2:end);
   dQ_lookup = Q_lookup(1:end - 1) - Q_lookup(2:end);

   % Collapse to the thermal-grid source term to show the quantity that the
   % main enthalpy solver actually consumes.
   [Sc_functions, chi_functions] = SPECTRALSOURCETERM(s.swd, s.albedo, ...
      s.I0, s.dz_spect, s.tau_N, s.tau_S, s.solar_dwavel, s.dz, ro_sno, ...
      s.z_nodes, s.z_nodes_spect, s.z_edges_spect, k_lookup_empty);

   [Sc_lookup, chi_lookup] = SPECTRALSOURCETERM(s.swd, s.albedo, ...
      s.I0, s.dz_spect, s.tau_N, s.tau_S, s.solar_dwavel, s.dz, ro_sno, ...
      s.z_nodes, s.z_nodes_spect, s.z_edges_spect, k_lookup);

   profiles = struct();
   profiles.z_nodes_spect = s.z_nodes_spect;
   profiles.z_flux_nodes = s.z_edges_spect;
   profiles.z_nodes_therm = s.z_nodes;
   profiles.bulk_functions = bulk_functions(1:numel(s.z_nodes_spect));
   profiles.bulk_lookup = bulk_lookup(1:numel(s.z_nodes_spect));
   profiles.Q_functions = Q_functions;
   profiles.Q_lookup = Q_lookup;
   profiles.dQ_functions = dQ_functions;
   profiles.dQ_lookup = dQ_lookup;
   profiles.Sc_functions = Sc_functions(:);
   profiles.Sc_lookup = Sc_lookup(:);
   profiles.chi_functions = chi_functions;
   profiles.chi_lookup = chi_lookup;
   profiles.up_functions = up_functions;
   profiles.dn_functions = dn_functions;
   profiles.up_lookup = up_lookup;
   profiles.dn_lookup = dn_lookup;
end

function exportTwoPanelPlot(z, exact, lookup, outfile, ttl, xlbl, ylbl, kwargs)
   %EXPORTTWOPANELPLOT Plot functions and lookup profiles with a difference panel.

   arguments
      z
      exact
      lookup
      outfile
      ttl
      xlbl
      ylbl
      kwargs.use_log_depth (1, 1) logical = false
   end

   % Plot the exact and lookup profiles together so the top panel shows the
   % shape agreement directly.
   f = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 760 820]);
   tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
   z_plot = z;
   if kwargs.use_log_depth
      z_plot = max(z_plot, min(z_plot(z_plot > 0)) / 10);
   end

   nexttile
   plot(exact, z_plot, 'k-', 'LineWidth', 1.8)
   hold on
   plot(lookup, z_plot, 'r--', 'LineWidth', 1.4)
   set(gca, 'YDir', 'reverse')
   if kwargs.use_log_depth
      set(gca, 'YScale', 'log')
   end
   xlabel(xlbl)
   ylabel(ylbl)
   title(ttl)
   legend({'functions', 'lookup'}, 'Location', 'best')
   grid on

   % Plot the lookup-minus-functions error separately so small differences
   % remain visible even when the top curves overlap.
   nexttile
   plot(lookup - exact, z_plot, 'b-', 'LineWidth', 1.4)
   set(gca, 'YDir', 'reverse')
   if kwargs.use_log_depth
      set(gca, 'YScale', 'log')
   end
   xlabel('lookup - functions')
   ylabel(ylbl)
   title('Difference')
   grid on

   exportgraphics(f, char(outfile), 'Resolution', 180)
   close(f)
end
