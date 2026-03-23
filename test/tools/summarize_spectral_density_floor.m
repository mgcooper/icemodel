function report = summarize_spectral_density_floor(kwargs)
   %SUMMARIZE_SPECTRAL_DENSITY_FLOOR Quantify usage of the ro>=300 clamp.
   %
   %  report = summarize_spectral_density_floor()
   %  report = summarize_spectral_density_floor(simyear=2016, solver=2)
   %  report = summarize_spectral_density_floor(output_file="/tmp/floor.mat")
   %
   % This diagnostic runs one formal icemodel case, reconstructs the
   % thermal-to-spectral density remap at every retained timestep, and reports:
   %  1. how often the spectral density floor is active
   %  2. the minimum raw remapped spectral density
   %  3. the worst bulk-extinction difference between floored and unfloored
   %     densities on those retained states

   arguments
      kwargs.simyear (1, 1) double {mustBeInteger, mustBePositive} = 2016
      kwargs.solver (1, 1) double {mustBeMember(kwargs.solver, [1 2 3])} = 2
      kwargs.smoke_site (1, :) string = "kanm"
      kwargs.output_file (1, :) string = ""
   end

   % Install the canonical suite config once so the formal case uses the same
   % environment as the accepted regression/perf tooling.
   [~, ~, ~, ~, suite_cleanup] = ...
      icemodel.test.helpers.bootstrapTestEnvironment(); %#ok<ASGLU>

   % Build and run one canonical formal smoke case, retaining the full formal
   % two-year contract so the density-floor audit sees the spinup evolution.
   cases = icemodel.test.helpers.getPerfCaseMatrix( ...
      tier="smoke", smbmodel="icemodel", solver=kwargs.solver, ...
      simyear=kwargs.simyear, smoke_sites=kwargs.smoke_site);
   if height(cases) ~= 1
      error('spectral density-floor summary expected one smoke case')
   end
   opts = icemodel.test.helpers.setModelOptsForCase(cases(1, :));
   opts = icemodel.resetopts(opts, 'output_profile', 'standard');
   [~, ice2, opts] = icemodel.test.helpers.runSmbModel(opts);

   % Rebuild the spectral grid once so each timestep reuses the same geometry.
   [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   [~, dz_spect, z_nodes_spect, ~, tau_N, tau_S, solar_dwavel] ...
      = EXTCOEFSINIT(opts, ro_ice);
   [~, ~, z_nodes_therm, ~] = CVMESH(opts.z0_thermal, opts.dz_thermal);

   % Scan all retained timesteps so the report quantifies both how often the
   % floor triggers and how large the induced bulk-coefficient perturbation is.
   min_raw_density = inf;
   total_nodes = 0;
   n_nodes_below_floor = 0;
   n_steps_with_floor = 0;
   worst_bulk_rel = 0;
   worst_step = 0;

   for istep = 1:size(ice2.f_ice, 2)
      ro_sno = ro_ice * ice2.f_ice(:, istep) + ro_liq * ice2.f_liq(:, istep);
      raw = interp1(z_nodes_therm, ro_sno, z_nodes_spect, ...
         'nearest', 'extrap');
      min_raw_density = min(min_raw_density, min(raw));

      below_floor = raw < 300;
      total_nodes = total_nodes + numel(raw);
      n_nodes_below_floor = n_nodes_below_floor + nnz(below_floor);

      if ~any(below_floor)
         continue
      end
      n_steps_with_floor = n_steps_with_floor + 1;

      % Compute k_bulk w and w/o the density floor
      k_bulk_floor = BULKEXTCOEFS(dz_spect, max(raw, 300), tau_N, tau_S, ...
         solar_dwavel);
      k_bulk_raw = computeBulkExactNoFloor(dz_spect, raw, tau_N, tau_S, ...
         solar_dwavel);

      rel = max(abs(k_bulk_floor - k_bulk_raw) ./ max(abs(k_bulk_floor), 1e-12));
      if rel > worst_bulk_rel
         worst_bulk_rel = rel;
         worst_step = istep;
      end
   end

   report = struct();
   report.case_id = string(cases.case_id(1));
   report.simyears = opts.simyears;
   report.n_spinup_years = opts.n_spinup_years;
   report.min_raw_density = min_raw_density;
   report.pct_nodes_below_floor = 100 * n_nodes_below_floor / total_nodes;
   report.pct_steps_with_floor = 100 * n_steps_with_floor / size(ice2.f_ice, 2);
   report.worst_bulk_rel = worst_bulk_rel;
   report.worst_step = worst_step;
   report.timestamp_utc = datetime('now', 'TimeZone', 'UTC');

   if ~isblanktext(kwargs.output_file)
      outdir = fileparts(char(kwargs.output_file));
      if ~isempty(outdir) && exist(outdir, 'dir') ~= 7
         mkdir(outdir);
      end
      save(char(kwargs.output_file), 'report');
   end

   report = struct2table(report);

   disp('Spectral density-floor summary:')
   disp(report(:, ...
      {'case_id', 'simyears', 'min_raw_density'}))
   disp(' ')
   disp(report(:, ...
      {'pct_nodes_below_floor', 'pct_steps_with_floor', 'worst_bulk_rel'}))
end

function bulkcoefs = computeBulkExactNoFloor(dz_spect, ro_sno, ...
      tau_N, tau_S, solar_dwavel)
   %COMPUTEBULKEXACTNOFLOOR Reproduce the exact transform without ro flooring.

   bulkcoefs = -log((sum(solar_dwavel .* exp(tau_S .* ro_sno), 2)) ...
      ./ (sum(solar_dwavel .* exp(tau_N .* ro_sno), 2))) / dz_spect(1);
   bulkcoefs = [bulkcoefs; bulkcoefs(end); bulkcoefs(end)];
end
