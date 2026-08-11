function tests = test_met_contracts
   %TEST_MET_CONTRACTS Verify met loading, processing, and saved-result paths.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   % Build one synthetic workspace and reuse it across the met contract
   % tests so the fixture cost stays out of each individual assertion.

   testCase.TestData.workspace = icemodel.test.fixtures.makeSyntheticWorkspace( ...
      [2015; 2016], configure=true, nsteps=24, dt_seconds=3600);
end

function teardown(testCase)
   % Remove the shared synthetic workspace after the file-level tests end.

   icemodel.test.fixtures.cleanupSyntheticWorkspace( ...
      testCase.TestData.workspace);
end

function test_loadmet_concatenates_years_and_computes_exchange(testCase)
   % Multi-year loads should concatenate years and derive exchange terms
   % like De exactly once through the namespaced loading path.

   workspace = testCase.TestData.workspace;
   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', [2015 2016]);

   [met, opts_out] = icemodel.loadmet(opts);

   testCase.verifyEqual(unique(year(met.Time))', [2015 2016]);
   testCase.verifyEqual(height(met), 48);
   testCase.verifyEqual(numel(opts_out.metfname), 2);
   testCase.verifyTrue(isvariable('De', met));
   testCase.verifyTrue(all(isfinite(met.De)));
   testCase.verifyGreaterThan(min(met.De), 0);

   [met_2015, opts_2015] = icemodel.loadmet(opts, 1);
   testCase.verifyEqual(unique(year(met_2015.Time)).', 2015);
   testCase.verifyEqual(opts_2015.simyears, 2015);
end

function test_loadmet_gates_promice_filled_by_window_coverage(testCase)
   % The canonical filled product cannot enter a run unless the producer
   % manifest pins every runtime artifact and the filled files' samples
   % cover every requested timestep with the required channels (POLICY
   % A4). Calendar-year ledger verdicts are producer bookkeeping and never
   % gate the runtime; geometry follows the POLICY A3 fallback chain and
   % never blocks loading.
   workspace = icemodel.test.fixtures.makeSyntheticWorkspace( ...
      [2015; 2016], configure=false, nsteps=96, dt_seconds=900);
   cleaner = onCleanup(@() ...
      icemodel.test.fixtures.cleanupSyntheticWorkspace(workspace));
   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', [2015 2016]);
   opts.forcings = 'promice_filled';
   opts.userdata = 'promice_filled';
   opts.readiness_file = fullfile(workspace.metdir, 'readiness.csv');
   opts.report_inputs_file = fullfile(workspace.metdir, ...
      'report-inputs.json');

   % Without the producer ledger artifact the gate refuses immediately.
   testCase.verifyError(@() icemodel.loadmet(opts), ...
      'icemodel:loadmet:promiceFilledNotReady');

   % Write a producer-schema ledger whose verdicts are all non-ready plus
   % the pinning manifest. Verdict strings do not gate the runtime: the
   % dateless request spans two full calendar years that the one-day
   % fixtures cannot cover, so the refusal is the window-coverage gate.
   site = repmat(string(opts.sitename), 2, 1);
   verdict_icemodel = ["not_forcing_ready"; "not_forcing_ready"];
   verdict_snowmodel = ["not_forcing_ready"; "not_forcing_ready"];
   ledger = table(site, [2015; 2016], verdict_icemodel, ...
      verdict_snowmodel, 'VariableNames', ...
      {'site', 'year', 'verdict_icemodel', 'verdict_snowmodel'});
   writetable(ledger, opts.readiness_file);
   writePromiceRuntimeManifest(opts);
   testCase.verifyError(@() icemodel.loadmet(opts), ...
      'icemodel:loadmet:promiceFilledWindowUncovered');

   % Narrowing the request to the covered fixture day passes the coverage
   % gate (arbitrary windows are first-class, POLICY A4) and moves the
   % refusal to product identity: these fixtures are not filled artifacts.
   opts.simyears = 2015;
   opts.numyears = 1;
   opts.startdate = datetime(2015, 1, 1, 0, 0, 0, 'TimeZone', 'UTC');
   opts.enddate = datetime(2015, 1, 1, 23, 0, 0, 'TimeZone', 'UTC');
   testCase.verifyError(@() icemodel.loadmet(opts), ...
      'icemodel:loadmet:promiceFilledIdentityMismatch');

   % Re-save the synthetic payloads under the exact filled-product
   % identity. The two files become CONTIGUOUS days of 2015 so one
   % requested window is fully covered across the file seam, matching the
   % contiguous span real producers publish (POLICY A6).
   canonical_files = cell(size(opts.metfname));
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   defaults = icemodel.forcing.reconstruct.setopts();
   provenance_channels = unique([defaults.plan_channels, ...
      icemodel.forcing.helpers.precipitationVariables(), ...
      "boom_height"], 'stable');
   for k = 1:numel(opts.metfname)
      S = load(opts.metfname{k}, 'met');
      met = S.met;
      % Day k of January 2015 at the fixture cadence replaces the
      % per-year fixture times so the two files abut seamlessly.
      met.Time = datetime(2015, 1, k, 0, 0, 0, 'TimeZone', 'UTC') ...
         + seconds(workspace.dt_seconds) * (0:height(met) - 1)';
      met.swu = met.albedo .* met.swd;
      met.rainf = zeros(height(met), 1);
      met.snowf = zeros(height(met), 1);
      met.boom_height = nan(height(met), 1);
      for channel = provenance_channels
         provenance = repmat(codes.observed, height(met), 1);
         provenance(~isfinite(met.(channel))) = codes.missing;
         met.(channel + "_provenance") = provenance;
      end
      metadata = met.Properties.UserData;
      metadata.site = string(opts.sitename);
      metadata.gapfill_product = "promice_filled";
      metadata.gapfill_engine_version = string(icemodel.internal.version());
      metadata.gapfill_policy_sha256 = ...
         icemodel.forcing.reconstruct.policySha256();
      metadata.gapfill_donors = string.empty(1, 0);
      metadata.gapfill_channels = defaults.plan_channels;
      metadata.gapfill_registry = codes;
      met.Properties.UserData = metadata;
      canonical_files{k} = fullfile(workspace.metdir, sprintf( ...
         'met_%s_promice_filled_201501%02d_201501%02d_15m.mat', ...
         opts.sitename, k, k));
      save(canonical_files{k}, 'met');
   end
   opts.metfname = canonical_files;
   % The request now spans both contiguous fixture days.
   opts.enddate = datetime(2015, 1, 2, 23, 45, 0, 'TimeZone', 'UTC');
   writePromiceRuntimeManifest(opts);
   S = load(opts.metfname{1}, 'met');
   canonical_met = S.met;
   met = canonical_met;
   met.Properties.UserData.gapfill_policy_sha256 = ...
      string(repmat('a', 1, 64));
   save(opts.metfname{1}, 'met');
   writePromiceRuntimeManifest(opts);
   testCase.verifyError(@() icemodel.loadmet(opts), ...
      'icemodel:loadmet:promiceFilledIdentityMismatch');
   met = canonical_met;
   met.Properties.UserData = rmfield( ...
      met.Properties.UserData, 'gapfill_policy_sha256');
   save(opts.metfname{1}, 'met');
   writePromiceRuntimeManifest(opts);
   testCase.verifyError(@() icemodel.loadmet(opts), ...
      'icemodel:loadmet:promiceFilledIdentityMismatch');
   met = canonical_met;
   met.Properties.UserData.gapfill_product = "promice";
   save(opts.metfname{1}, 'met');
   writePromiceRuntimeManifest(opts);
   testCase.verifyError(@() icemodel.loadmet(opts), ...
      'icemodel:loadmet:promiceFilledIdentityMismatch');
   met = canonical_met;
   met.Properties.UserData.site = "other";
   save(opts.metfname{1}, 'met');
   writePromiceRuntimeManifest(opts);
   testCase.verifyError(@() icemodel.loadmet(opts), ...
      'icemodel:loadmet:promiceFilledIdentityMismatch');
   met = canonical_met;
   save(opts.metfname{1}, 'met');
   writePromiceRuntimeManifest(opts);
   % A forged _15m filename cannot make an hourly payload runnable.
   met = canonical_met(1:4:end, :);
   save(opts.metfname{1}, 'met');
   writePromiceRuntimeManifest(opts);
   testCase.verifyError(@() icemodel.loadmet(opts), ...
      'icemodel:loadmet:promiceFilledCadenceMismatch');
   met = canonical_met;
   save(opts.metfname{1}, 'met');
   writePromiceRuntimeManifest(opts);
   met = canonical_met;
   met.Properties.UserData = rmfield( ...
      met.Properties.UserData, 'gapfill_registry');
   save(opts.metfname{1}, 'met');
   writePromiceRuntimeManifest(opts);
   testCase.verifyError(@() icemodel.loadmet(opts), ...
      'icemodel:loadmet:promiceFilledProvenanceMismatch');
   met = canonical_met;
   met.tair_provenance = [];
   save(opts.metfname{1}, 'met');
   writePromiceRuntimeManifest(opts);
   testCase.verifyError(@() icemodel.loadmet(opts), ...
      'icemodel:loadmet:promiceFilledProvenanceMismatch');
   met = canonical_met;
   save(opts.metfname{1}, 'met');
   writePromiceRuntimeManifest(opts);

   % A partial source split remains honest only when every finite phase can
   % fit inside the finite total. The snow runtime gate must reject an
   % impossible component instead of accepting the other valid component.
   snow_opts = opts;
   snow_opts.smbmodel = 'snowmodel';
   met = canonical_met;
   met.rainf(1) = met.ppt(1) + 1;
   met.snowf(1) = NaN;
   met.snowf_provenance(1) = codes.missing;
   save(opts.metfname{1}, 'met');
   writePromiceRuntimeManifest(opts);
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.verifyPromiceFilledReadiness(snow_opts, []), ...
      'icemodel:loadmet:promiceFilledWindowUncovered');
   met = canonical_met;
   met.rainf(1) = NaN;
   met.rainf_provenance(1) = codes.missing;
   met.snowf(1) = met.ppt(1) + 1;
   save(opts.metfname{1}, 'met');
   writePromiceRuntimeManifest(opts);
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.verifyPromiceFilledReadiness(snow_opts, []), ...
      'icemodel:loadmet:promiceFilledWindowUncovered');
   % A nonnegative partial split at or below total remains admissible.
   met = canonical_met;
   met.rainf(1) = 0.5 * met.ppt(1);
   met.snowf(1) = NaN;
   met.snowf_provenance(1) = codes.missing;
   save(opts.metfname{1}, 'met');
   writePromiceRuntimeManifest(opts);
   testCase.verifyWarningFree(@() ...
      icemodel.forcing.reconstruct.verifyPromiceFilledReadiness(snow_opts, []));
   met = canonical_met;
   save(opts.metfname{1}, 'met');
   writePromiceRuntimeManifest(opts);

   % Geometry never blocks a run or readiness (POLICY A3 / D-0a): with no
   % valid measured boom sample anywhere, loading succeeds on the nominal
   % single-source constant with one warning, and the resolution is
   % recorded in the returned options.
   nominal = icemodel.parameterLookup('promice_nominal_boom_height_m');
   [met_nominal, opts_nominal] = testCase.verifyWarning( ...
      @() icemodel.loadmet(opts), ...
      'icemodel:loadmet:promiceBoomHeightNominal');
   testCase.verifyEqual(opts_nominal.z_tair, ...
      repmat(nominal, height(met_nominal), 1));
   testCase.verifyEqual(opts_nominal.boom_height_source, 'nominal');
   testCase.verifyEqual(opts_nominal.boom_height_fraction_fallback, 1);

   % Boom height is optional product geometry, so an artifact may omit both
   % its value and provenance columns and still use the same nominal runtime
   % fallback as an all-missing column.
   for k = 1:numel(opts.metfname)
      filename = opts.metfname{k};
      S = load(filename, 'met');
      met = removevars(S.met, {'boom_height', 'boom_height_provenance'});
      save(filename, 'met');
   end
   writePromiceRuntimeManifest(opts);
   [met_nominal, opts_nominal] = testCase.verifyWarning( ...
      @() icemodel.loadmet(opts), ...
      'icemodel:loadmet:promiceBoomHeightNominal');
   testCase.verifyEqual(opts_nominal.z_tair, ...
      repmat(nominal, height(met_nominal), 1));
   testCase.verifyEqual(opts_nominal.boom_height_source, 'nominal');
   testCase.verifyEqual(opts_nominal.boom_height_fraction_fallback, 1);

   % Restore the optional columns for the measured-geometry checks below.
   for k = 1:numel(opts.metfname)
      filename = opts.metfname{k};
      S = load(filename, 'met');
      met = S.met;
      met.boom_height = nan(height(met), 1);
      met.boom_height_provenance = ...
         repmat(codes.missing, height(met), 1);
      save(filename, 'met');
   end
   writePromiceRuntimeManifest(opts);

   % Install one changing boom-height series in each source file. The load
   % contract returns it as all three observation heights and computes De from
   % that series rather than the unresolved setopts sentinel.
   expected_parts = cell(numel(opts.metfname), 1);
   for k = 1:numel(opts.metfname)
      filename = opts.metfname{k};
      S = load(filename, 'met');
      met = S.met;
      met.boom_height = linspace(2.6, 2.2, height(met)).';
      met.boom_height_provenance(:) = codes.observed;
      expected_parts{k} = met.boom_height;
      save(filename, 'met');
   end
   expected = vertcat(expected_parts{:});
   writePromiceRuntimeManifest(opts);
   [met, opts_out] = icemodel.loadmet(opts);
   expected_height = height(met);
   testCase.verifyEqual(unique(year(met.Time)).', 2015);
   % A fully measured series is recorded with a zero fallback fraction.
   testCase.verifyEqual(opts_out.boom_height_source, 'measured');
   testCase.verifyEqual(opts_out.boom_height_fraction_fallback, 0);
   testCase.verifyTrue(opts_out.promice_filled_readiness_verified);
   testCase.verifyTrue(opts_out.promice_filled_manifest_verified);
   testCase.verifyEqual(opts_out.promice_filled_verified_forcing, ...
      opts_out.forcings);
   testCase.verifyEqual(opts_out.promice_filled_verified_site, ...
      lower(opts_out.sitename));
   testCase.verifyEqual(opts_out.promice_filled_verified_simyears(:), ...
      opts_out.simyears(:));
   testCase.verifyEqual(opts_out.promice_filled_verified_metfname(:), ...
      opts_out.metfname(:));
   testCase.verifyEqual(opts_out.promice_filled_verified_startdate, ...
      opts_out.startdate);
   testCase.verifyEqual(opts_out.promice_filled_verified_enddate, ...
      opts_out.enddate);
   testCase.verifyEqual(opts_out.promice_filled_verified_calendar_type, ...
      opts_out.calendar_type);
   testCase.verifyEqual(opts_out.promice_filled_verified_smbmodel, ...
      opts_out.smbmodel);
   testCase.verifyEqual(opts_out.promice_filled_verified_dt, opts_out.dt);
   testCase.verifyTrue( ...
      icemodel.forcing.reconstruct.promiceFilledVerificationMatches(opts_out));
   % A requested file subset verifies only the bytes it will load. An
   % unselected file whose bytes do not match the producer manifest
   % cannot invalidate a valid partial load.
   full_met = met;
   second_payload = load(opts.metfname{2}, 'met');
   second_intact = second_payload.met;
   met = second_intact;
   met.tair(1) = met.tair(1) + 1;
   save(opts.metfname{2}, 'met');
   first_payload = load(opts.metfname{1}, 'met');
   [partial_met, partial_out] = icemodel.loadmet(opts, 1);
   testCase.verifyEqual(partial_met.Time, first_payload.met.Time);
   testCase.verifyEqual(partial_out.startdate, min(first_payload.met.Time));
   testCase.verifyEqual(partial_out.enddate, max(first_payload.met.Time));
   testCase.verifyTrue( ...
      icemodel.forcing.reconstruct.promiceFilledVerificationMatches(partial_out, 1));
   testCase.verifyFalse( ...
      icemodel.forcing.reconstruct.promiceFilledVerificationMatches(partial_out));
   met = second_intact;
   save(opts.metfname{2}, 'met');
   met = full_met;
   changed = opts_out;
   changed.promice_filled_readiness_verified = false;
   testCase.verifyFalse( ...
      icemodel.forcing.reconstruct.promiceFilledVerificationMatches(changed));
   changed = opts_out;
   changed.promice_filled_manifest_verified = false;
   testCase.verifyFalse( ...
      icemodel.forcing.reconstruct.promiceFilledVerificationMatches(changed));
   changed = opts_out;
   changed.sitename = 'other';
   testCase.verifyFalse( ...
      icemodel.forcing.reconstruct.promiceFilledVerificationMatches(changed));
   changed = opts_out;
   changed.simyears = [changed.simyears(:).', 2016];
   testCase.verifyFalse( ...
      icemodel.forcing.reconstruct.promiceFilledVerificationMatches(changed));
   changed = opts_out;
   changed.metfname{1} = [changed.metfname{1} '.other'];
   testCase.verifyFalse( ...
      icemodel.forcing.reconstruct.promiceFilledVerificationMatches(changed));
   changed = opts_out;
   changed.forcings = 'promice';
   testCase.verifyFalse( ...
      icemodel.forcing.reconstruct.promiceFilledVerificationMatches(changed));
   testCase.verifyError(@() icemodel.loadmet(changed), ...
      'icemodel:loadmet:promiceFilledIdentityMismatch');
   changed = opts_out;
   changed.enddate = changed.enddate - minutes(15);
   testCase.verifyFalse( ...
      icemodel.forcing.reconstruct.promiceFilledVerificationMatches(changed));
   changed = opts_out;
   changed.calendar_type = 'standard';
   testCase.verifyFalse( ...
      icemodel.forcing.reconstruct.promiceFilledVerificationMatches(changed));
   changed = opts_out;
   changed.smbmodel = 'snowmodel';
   testCase.verifyFalse( ...
      icemodel.forcing.reconstruct.promiceFilledVerificationMatches(changed));
   changed = opts_out;
   changed.dt = 1800;
   testCase.verifyFalse( ...
      icemodel.forcing.reconstruct.promiceFilledVerificationMatches(changed));
   testCase.verifyError(@() icemodel.loadmet(changed), ...
      'icemodel:loadmet:promiceFilledCadenceMismatch');
   testCase.verifyEqual(opts_out.z_tair, expected, 'AbsTol', 1e-12);
   testCase.verifyEqual(opts_out.z_wind, expected, 'AbsTol', 1e-12);
   testCase.verifyEqual(opts_out.z_relh, expected, 'AbsTol', 1e-12);
   expected_De = ...
      icemodel.surface.turbulence.bulk_richardson.exchange_coefficients( ...
      met.wspd, opts_out.z0_bulk, expected, expected);
   testCase.verifyEqual(met.De, expected_De, 'AbsTol', 1e-12);
   % A custom producer plan cannot omit typed provenance for a channel
   % the runtime still requires (POLICY A4/A7).
   filename = opts.metfname{1};
   S = load(filename, 'met');
   intact = S.met;
   met = intact;
   met.Properties.UserData.gapfill_channels = "tair";
   met = removevars(met, "lwd_provenance");
   save(filename, 'met');
   writePromiceRuntimeManifest(opts);
   testCase.verifyError(@() icemodel.loadmet(opts), ...
      'icemodel:loadmet:promiceFilledProvenanceMismatch');
   met = intact;
   save(filename, 'met');
   writePromiceRuntimeManifest(opts);
   backup = fullfile(workspace.metdir, 'filled-backup.mat');
   copyfile(filename, backup);
   S = load(filename, 'met');
   tampered = S.met;
   tampered.tair(1) = tampered.tair(1) + 1;
   met = tampered;
   save(filename, 'met');
   testCase.verifyError(@() icemodel.loadmet(opts), ...
      'icemodel:loadmet:promiceFilledIdentityMismatch');
   copyfile(backup, filename, 'f');
   [ice1, ~, run_opts] = skinmodel(opts);
   testCase.verifyEqual(numel(ice1.Tsfc), expected_height);
   testCase.verifyEqual(run_opts.z_tair, expected, 'AbsTol', 1e-12);
   % Legacy labels cannot be used to bypass the filled-product runtime gate.
   for alias = ["kanm", "kanl"]
      alias_opts = opts;
      alias_opts.forcings = alias;
      alias_opts.userdata = alias;
      testCase.verifyError(@() icemodel.loadmet(alias_opts), ...
         'icemodel:loadmet:promiceFilledIdentityMismatch');
   end

   % A nonfinite sample and a measured value at or below z0_bulk do not
   % block the run even when the ledger is stale: both demote into the
   % interpolated rung with one warning (POLICY A3), and the bridged
   % series feeds every observation height.
   filename = opts.metfname{1};
   S = load(filename, 'met');
   met_invalid = S.met;
   met_invalid.boom_height(2) = NaN;
   met_invalid.boom_height_provenance(2) = codes.missing;
   met_invalid.boom_height(4) = opts.z0_bulk / 2;
   met = met_invalid;
   save(filename, 'met');
   writePromiceRuntimeManifest(opts);
   [~, opts_interp] = testCase.verifyWarning(@() icemodel.loadmet(opts), ...
      'icemodel:loadmet:promiceBoomHeightInterpolated');
   testCase.verifyEqual(opts_interp.boom_height_source, 'interpolated');
   testCase.verifyEqual(opts_interp.boom_height_fraction_fallback, ...
      2 / expected_height);
   % Interior invalid samples bridge linearly between valid neighbors.
   expected_bridge = expected;
   expected_bridge(2) = (expected(1) + expected(3)) / 2;
   expected_bridge(4) = (expected(3) + expected(5)) / 2;
   testCase.verifyEqual(opts_interp.z_tair, expected_bridge, ...
      'AbsTol', 1e-12);
   alias_opts.forcings = "kanm";
   alias_opts.userdata = "kanm";
   testCase.verifyError(@() icemodel.loadmet(alias_opts), ...
      'icemodel:loadmet:promiceFilledIdentityMismatch');

   % Restore intact geometry, then punch one required-channel hole inside
   % the window: the gate refuses with the named coverage error naming the
   % channel compactly, and narrowing the window past the hole loads again.
   met = met_invalid;
   met.boom_height = expected(1:height(met_invalid));
   met.boom_height_provenance(:) = codes.observed;
   save(filename, 'met');
   hole_file = opts.metfname{2};
   S = load(hole_file, 'met');
   met_intact2 = S.met;
   met_hole = met_intact2;
   met_hole.tair(38) = NaN;
   met_hole.tair_provenance(38) = codes.missing;
   met = met_hole;
   save(hole_file, 'met');
   writePromiceRuntimeManifest(opts);
   try
      icemodel.loadmet(opts);
      testCase.verifyFail('expected a window-coverage refusal');
   catch ME
      testCase.verifyEqual(ME.identifier, ...
         'icemodel:loadmet:promiceFilledWindowUncovered');
      % The refusal names the uncovered channel compactly.
      testCase.verifyTrue(contains(ME.message, 'tair'));
   end

   % A partial-window request that stops before the hole is first-class
   % (POLICY A4 water years / arbitrary windows) and loads: coverage is
   % judged against the requested window only.
   narrowed = opts;
   narrowed.enddate = datetime(2015, 1, 2, 8, 0, 0, 'TimeZone', 'UTC');
   met_window = icemodel.loadmet(narrowed);
   testCase.verifyEqual(height(met_window), 129);
   testCase.verifyEqual(max(met_window.Time), narrowed.enddate);

   % A physically absent timestep row violates the published contiguous
   % 15-minute cadence even when every present sample is finite.
   met = met_intact2;
   met(10, :) = [];
   save(hole_file, 'met');
   writePromiceRuntimeManifest(opts);
   try
      icemodel.loadmet(opts);
      testCase.verifyFail('expected a missing-timestep refusal');
   catch ME
      testCase.verifyEqual(ME.identifier, ...
         'icemodel:loadmet:promiceFilledCadenceMismatch');
      testCase.verifyTrue(contains(ME.message, 'contiguous 15-minute'));
   end

   % Removing the first row leaves a contiguous and correctly phased file,
   % but cannot redefine the requested grid to start at 00:15.
   met = met_intact2;
   save(hole_file, 'met');
   first_file = opts.metfname{1};
   S = load(first_file, 'met');
   met_intact1 = S.met;
   met = met_intact1;
   met(1, :) = [];
   save(first_file, 'met');
   writePromiceRuntimeManifest(opts);
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.verifyPromiceFilledReadiness(opts), ...
      'icemodel:loadmet:promiceFilledWindowUncovered');

   % A complete but phase-shifted file is not the canonical 15-minute
   % product, even though its adjacent rows remain 15 minutes apart.
   met = met_intact1;
   met.Time = met.Time + minutes(5);
   save(first_file, 'met');
   writePromiceRuntimeManifest(opts);
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.verifyPromiceFilledReadiness(opts), ...
      'icemodel:loadmet:promiceFilledCadenceMismatch');
   met = met_intact1;
   save(first_file, 'met');

   % A finite but out-of-bounds required value is not runtime coverage:
   % readiness applies the A15 scalar registry while leaving relational
   % shortwave checks as producer-ledger diagnostics.
   met = met_intact2;
   met.lwd(10) = 500;
   save(hole_file, 'met');
   writePromiceRuntimeManifest(opts);
   try
      icemodel.forcing.reconstruct.verifyPromiceFilledReadiness(opts);
      testCase.verifyFail('expected a scalar-bound refusal');
   catch ME
      testCase.verifyEqual(ME.identifier, ...
         'icemodel:loadmet:promiceFilledWindowUncovered');
      testCase.verifyTrue(contains(ME.message, 'lwd'));
   end

   % Restore the intact second file, then prove the A5 snow-model
   % extension through the public verifier: a snowfall-consuming smbmodel
   % additionally requires finite total ppt OR snowf at every requested
   % timestep, while the historical zero-rain models gate only on the
   % seven-channel icemodel set.
   met = met_intact2;
   save(hole_file, 'met');
   writePromiceRuntimeManifest(opts);
   snow_opts = opts;
   snow_opts.smbmodel = 'snowmodel';
   returned = icemodel.forcing.reconstruct.verifyPromiceFilledReadiness(snow_opts);
   testCase.verifyTrue(returned.promice_filled_readiness_verified);
   % Scalar-valid precipitation is still unusable when a complete source
   % split violates A10; the runtime gate must fail before the default
   % source-phase resolver sees it.
   met_bad = met_intact2;
   met_bad.ppt(10) = 1e-6;
   met_bad.rainf(10) = 0.75e-6;
   met_bad.snowf(10) = 0.75e-6;
   met = met_bad;
   save(hole_file, 'met');
   writePromiceRuntimeManifest(opts);
   try
      icemodel.forcing.reconstruct.verifyPromiceFilledReadiness(snow_opts);
      testCase.verifyFail('expected a precipitation-identity refusal');
   catch ME
      testCase.verifyEqual(ME.identifier, ...
         'icemodel:loadmet:promiceFilledWindowUncovered');
      testCase.verifyTrue(contains(ME.message, 'ppt|snowf'));
   end
   % Remove the precipitation coverage with honest provenance: the
   % snow-model gate refuses; the skinmodel gate still passes because
   % rainf alone never satisfies (and never gates) either verdict set.
   met_dry = met_intact2;
   met_dry.ppt(:) = NaN;
   met_dry.ppt_provenance(:) = codes.missing;
   met_dry.snowf(:) = NaN;
   met_dry.snowf_provenance(:) = codes.missing;
   met = met_dry;
   save(hole_file, 'met');
   writePromiceRuntimeManifest(opts);
   snow_opts = opts;
   snow_opts.smbmodel = 'snowmodel';
   testCase.verifyError( ...
      @() icemodel.forcing.reconstruct.verifyPromiceFilledReadiness(snow_opts), ...
      'icemodel:loadmet:promiceFilledWindowUncovered');
   returned = icemodel.forcing.reconstruct.verifyPromiceFilledReadiness(opts);
   testCase.verifyTrue(returned.promice_filled_readiness_verified);
   clear cleaner
end

function test_overlapping_promice_filled_files_are_rejected(testCase)
   % The MATLAB verifier and loader must reject duplicate interval starts
   % before code-generation trust flags can describe the selected artifacts.
   [workspace, opts, ~] = makeCodegenPromiceFixture();
   cleaner = onCleanup(@() ...
      icemodel.test.fixtures.cleanupSyntheticWorkspace(workspace));
   original_file = opts.metfname{1};
   overlap_file = fullfile(workspace.metdir, sprintf( ...
      'met_%s_promice_filled_20150101_20150101_overlap_15m.mat', ...
      opts.sitename));
   copyfile(original_file, overlap_file)
   opts.metfname = {original_file, overlap_file};
   writePromiceRuntimeManifest(opts)

   identifier = 'icemodel:loadmet:promiceFilledIntervalOverlap';
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.verifyPromiceFilledReadiness(opts), identifier);
   testCase.verifyFalse(opts.promice_filled_readiness_verified);
   testCase.verifyFalse(opts.promice_filled_manifest_verified);
   testCase.verifyFalse(opts.promice_filled_provenance_verified);
   testCase.verifyError(@() icemodel.loadmet(opts), identifier);
   clear cleaner
end

function test_file_subset_normalizes_selected_times_before_bounds(testCase)
   % A selected file's staged coordinate must enter the UTC comparison frame
   % before fileiter narrows the verified request and loadmet filters rows.
   [workspace, opts, canonical_met] = makeCodegenPromiceFixture();
   cleaner = onCleanup(@() ...
      icemodel.test.fixtures.cleanupSyntheticWorkspace(workspace));
   selected_file = opts.metfname{1};
   unselected_file = fullfile(workspace.metdir, sprintf( ...
      'met_%s_promice_filled_20150102_20150102_15m.mat', ...
      opts.sitename));
   copyfile(selected_file, unselected_file)
   opts.metfname = {selected_file, unselected_file};

   variants = ["unzoned", "non-utc"];
   for variant = variants
      met = canonical_met;
      if variant == "unzoned"
         met.Time.TimeZone = '';
      else
         met.Time.TimeZone = 'America/New_York';
      end
      artifact_metadata = met.Properties.UserData;
      save(selected_file, 'met', 'artifact_metadata')
      writePromiceRuntimeManifest(opts)

      [loaded, verified] = icemodel.loadmet(opts, 1);
      testCase.verifyEqual(string(loaded.Time.TimeZone), "UTC");
      testCase.verifyTrue(all(loaded.Time == canonical_met.Time));
      testCase.verifyEqual(verified.startdate, opts.startdate);
      testCase.verifyEqual(verified.enddate, opts.enddate);
      testCase.verifyTrue( ...
         icemodel.forcing.reconstruct.promiceFilledVerificationMatches(verified, 1));
   end
   clear cleaner
end

function test_readiness_rejects_unsafe_station_token(testCase)
   % The default ledger path must not accept separators or traversal text.
   opts = struct('forcings', "promice_filled", ...
      'sitename', "../escape", 'simyears', 2020);

   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.verifyPromiceFilledReadiness(opts), ...
      'icemodel:reconstruct:mustBeStationToken:invalidToken');
end

function test_codegen_flags_require_canonical_artifact_provenance(testCase)
   % A ready ledger and a freshly repinned manifest cannot restore
   % generated-code trust after the artifact policy or the required
   % provenance becomes invalid.
   [workspace, opts, canonical_met] = makeCodegenPromiceFixture();
   cleaner = onCleanup(@() ...
      icemodel.test.fixtures.cleanupSyntheticWorkspace(workspace));
   filename = opts.metfname{1};

   % The canonical MATLAB preparation path sets all three trust flags, and
   % the codegen-compatible check in generated loading accepts that exact opts.
   verified = icemodel.forcing.reconstruct.verifyPromiceFilledReadiness(opts);
   testCase.verifyTrue(verified.promice_filled_readiness_verified);
   testCase.verifyTrue(verified.promice_filled_manifest_verified);
   testCase.verifyTrue(verified.promice_filled_provenance_verified);
   testCase.verifyTrue(icemodel.forcing.reconstruct.promiceFilledVerificationMatches(verified));
   legacy = verified;
   legacy.promice_filled_provenance_verified = false;
   testCase.verifyFalse(icemodel.forcing.reconstruct.promiceFilledVerificationMatches(legacy));

   % Re-pin each deliberate artifact corruption so the refusal proves the
   % provenance gate, not a stale producer-manifest digest. The ready ledger
   % remains unchanged and therefore cannot rescue either payload.
   variants = ["stale-policy", "missing-channel"];
   identifiers = ["icemodel:loadmet:promiceFilledIdentityMismatch", ...
      "icemodel:loadmet:promiceFilledProvenanceMismatch"];
   for k = 1:numel(variants)
      met = canonical_met;
      if variants(k) == "stale-policy"
         met.Properties.UserData.gapfill_policy_sha256 = ...
            string(repmat('a', 1, 64));
      else
         met = removevars(met, 'lwd_provenance');
      end
      artifact_metadata = met.Properties.UserData;
      save(filename, 'met', 'artifact_metadata')
      writePromiceRuntimeManifest(opts)

      testCase.verifyError(@() ...
         icemodel.forcing.reconstruct.verifyPromiceFilledReadiness(opts), identifiers(k));
      testCase.verifyFalse(opts.promice_filled_readiness_verified);
      testCase.verifyFalse(opts.promice_filled_manifest_verified);
      testCase.verifyFalse(opts.promice_filled_provenance_verified);
      testCase.verifyFalse(icemodel.forcing.reconstruct.promiceFilledVerificationMatches(opts));
   end

   % Restoring and repinning the canonical bytes retains the supported path.
   met = canonical_met;
   artifact_metadata = met.Properties.UserData;
   save(filename, 'met', 'artifact_metadata')
   writePromiceRuntimeManifest(opts)
   restored = icemodel.forcing.reconstruct.verifyPromiceFilledReadiness(opts);
   testCase.verifyTrue(icemodel.forcing.reconstruct.promiceFilledVerificationMatches(restored));
   clear cleaner
end

function test_initforcings_keeps_finite_precip_out_of_solver(testCase)
   % Reconstructed rain remains a data product and must not activate unfinished
   % precipitation physics in existing model runs.
   workspace = testCase.TestData.workspace;
   rainf = linspace(0, 1e-6, workspace.nsteps).';
   snowf = 2e-6 * ones(workspace.nsteps, 1);
   ppt = rainf + snowf;
   icemodel.test.fixtures.makeSyntheticMetFile(2016, ...
      'sitename', workspace.sitename, ...
      'forcings', workspace.forcings, ...
      'nsteps', workspace.nsteps, ...
      'dt_seconds', workspace.dt_seconds, ...
      'ppt', ppt, ...
      'rainf', rainf, ...
      'snowf', snowf, ...
      'metdir', workspace.metdir);
   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016);

   [~, ~, ~, ~, ~, ~, ~, returned] = ...
      icemodel.surface.initialize_surface_forcings(opts);

   testCase.verifyEqual(returned, zeros(workspace.nsteps, 1));
end

function test_initforcings_keeps_missing_precip_out_of_solver(testCase)
   % An allowed all-missing precipitation placeholder preserves the
   % historical no-rain fallback and cannot inject NaN into surface physics.
   workspace = testCase.TestData.workspace;
   icemodel.test.fixtures.makeSyntheticMetFile(2016, ...
      'sitename', workspace.sitename, ...
      'forcings', workspace.forcings, ...
      'nsteps', workspace.nsteps, ...
      'dt_seconds', workspace.dt_seconds, ...
      'ppt', nan(workspace.nsteps, 1), ...
      'rainf', nan(workspace.nsteps, 1), ...
      'snowf', nan(workspace.nsteps, 1), ...
      'metdir', workspace.metdir);
   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016);

   [~, ~, ~, ~, ~, ~, ~, returned] = ...
      icemodel.surface.initialize_surface_forcings(opts);

   testCase.verifyEqual(returned, zeros(workspace.nsteps, 1));
end

function test_loadmet_boom_chain_for_alias_forcing(testCase)
   % A legacy-alias met file carrying the boom channel resolves measured
   % heights, demotes invalid samples (nonfinite or at/below z0_bulk) into
   % the interpolated rung with one warning, and lands on the nominal
   % constant when nothing measured exists — geometry never errors
   % (POLICY A3 / D-0a).
   workspace = testCase.TestData.workspace;
   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016);
   filename = opts.metfname{1};

   % 1) Fully valid measured series: heights adopt the measurements with
   % no warning and a zero fallback fraction.
   S = load(filename, 'met');
   met_measured = S.met;
   measured = linspace(2.7, 2.4, height(met_measured)).';
   met_measured.boom_height = measured;
   met = met_measured;
   save(filename, 'met');
   [~, opts_measured] = testCase.verifyWarningFree( ...
      @() icemodel.loadmet(opts));
   testCase.verifyEqual(opts_measured.z_tair, measured, 'AbsTol', 1e-12);
   testCase.verifyEqual(opts_measured.boom_height_source, 'measured');
   testCase.verifyEqual(opts_measured.boom_height_fraction_fallback, 0);

   % 2) Invalid samples demote into the interpolated rung: interior
   % samples bridge linearly and the series ends extend the nearest valid
   % measurement, with one warning identifying the degradation.
   met_invalid = met_measured;
   met_invalid.boom_height(1) = NaN;             % leading-edge extension
   met_invalid.boom_height(5) = NaN;             % interior linear bridge
   met_invalid.boom_height(9) = opts.z0_bulk;    % z0 demotion (at bound)
   met_invalid.boom_height(end) = -1;            % trailing-edge extension
   met = met_invalid;
   save(filename, 'met');
   [~, opts_interp] = testCase.verifyWarning(@() icemodel.loadmet(opts), ...
      'icemodel:loadmet:promiceBoomHeightInterpolated');
   expected = measured;
   expected(1) = measured(2);
   expected(5) = (measured(4) + measured(6)) / 2;
   expected(9) = (measured(8) + measured(10)) / 2;
   expected(end) = measured(end - 1);
   testCase.verifyEqual(opts_interp.z_tair, expected, 'AbsTol', 1e-12);
   testCase.verifyEqual(opts_interp.z_wind, expected, 'AbsTol', 1e-12);
   testCase.verifyEqual(opts_interp.z_relh, expected, 'AbsTol', 1e-12);
   testCase.verifyEqual(opts_interp.boom_height_source, 'interpolated');
   testCase.verifyEqual(opts_interp.boom_height_fraction_fallback, ...
      4 / height(met_invalid));

   % 3) No valid measurement anywhere: the nominal single-source constant
   % applies to every sample with one warning.
   met_none = met_measured;
   met_none.boom_height(:) = NaN;
   met = met_none;
   save(filename, 'met');
   [~, opts_nominal] = testCase.verifyWarning(@() icemodel.loadmet(opts), ...
      'icemodel:loadmet:promiceBoomHeightNominal');
   nominal = icemodel.parameterLookup('promice_nominal_boom_height_m');
   testCase.verifyEqual(opts_nominal.z_tair, ...
      repmat(nominal, height(met_none), 1));
   testCase.verifyEqual(opts_nominal.boom_height_source, 'nominal');
   testCase.verifyEqual(opts_nominal.boom_height_fraction_fallback, 1);
end

function test_loadmet_alias_without_boom_channel_keeps_nominal(testCase)
   % Historical alias fixtures without the boom channel keep their scalar
   % nominal setopts geometry with no warning (the documented alias
   % contract, not a degradation); the resolution is still recorded.
   workspace = testCase.TestData.workspace;
   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016);
   nominal = icemodel.parameterLookup('promice_nominal_boom_height_m');
   [~, opts_out] = testCase.verifyWarningFree(@() icemodel.loadmet(opts));
   testCase.verifyEqual(opts_out.z_tair, nominal);
   testCase.verifyEqual(opts_out.boom_height_source, 'nominal');
   testCase.verifyEqual(opts_out.boom_height_fraction_fallback, 1);
end

function test_initforcings_phase_source_default_exposes_product_split(testCase)
   % The default 'source' phase option exposes the met product's own
   % rainf/snowf split bit-identically while the solver rain forcing stays
   % zero (POLICY A10 / D-18; zero-rain per D-0b).
   workspace = testCase.TestData.workspace;
   rainf = linspace(0, 1e-6, workspace.nsteps).';
   snowf = 2e-6 * ones(workspace.nsteps, 1);
   ppt = rainf + snowf;
   icemodel.test.fixtures.makeSyntheticMetFile(2016, ...
      'sitename', workspace.sitename, ...
      'forcings', workspace.forcings, ...
      'nsteps', workspace.nsteps, ...
      'dt_seconds', workspace.dt_seconds, ...
      'ppt', ppt, ...
      'rainf', rainf, ...
      'snowf', snowf, ...
      'metdir', workspace.metdir);
   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016);

   [~, ~, ~, ~, ~, ~, ~, rain, ~, ~, ~, ~, returned_rainf, ...
      returned_snowf] = icemodel.surface.initialize_surface_forcings(opts);

   testCase.verifyEqual(rain, zeros(workspace.nsteps, 1));
   % Bit-identical to the shipped components: no tolerance.
   testCase.verifyEqual(returned_rainf, rainf);
   testCase.verifyEqual(returned_snowf, snowf);
end

function test_initforcings_phase_threshold_partitions_total(testCase)
   % The 'threshold' phase option repartitions the canonical total by air
   % temperature with the single-source transition temperature; the
   % shipped split is ignored in this mode and conservation holds by
   % construction (POLICY A10 / D-18).
   workspace = testCase.TestData.workspace;
   transition = icemodel.forcing.reconstruct.setopts(). ...
      rain_snow_transition_temperature_k;
   ppt = 1e-6 * ones(workspace.nsteps, 1);
   % A shipped all-rain split that disagrees with the threshold result
   % proves the option selects the partition, not the product components.
   icemodel.test.fixtures.makeSyntheticMetFile(2016, ...
      'sitename', workspace.sitename, ...
      'forcings', workspace.forcings, ...
      'nsteps', workspace.nsteps, ...
      'dt_seconds', workspace.dt_seconds, ...
      'ppt', ppt, ...
      'rainf', ppt, ...
      'snowf', zeros(workspace.nsteps, 1), ...
      'metdir', workspace.metdir);
   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016);
   % Straddle the transition so both phases appear in the partition.
   filename = opts.metfname{1};
   S = load(filename, 'met');
   met_straddle = S.met;
   met_straddle.tair(1:2:end) = transition + 1;
   met_straddle.tair(2:2:end) = transition - 1;
   met = met_straddle;
   save(filename, 'met');
   opts = icemodel.resetopts(opts, 'precip_phase_source', 'threshold');

   [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, returned_rainf, ...
      returned_snowf] = icemodel.surface.initialize_surface_forcings(opts);

   expected_rainf = ppt .* (met_straddle.tair >= transition);
   expected_snowf = ppt .* (met_straddle.tair < transition);
   testCase.verifyEqual(returned_rainf, expected_rainf);
   testCase.verifyEqual(returned_snowf, expected_snowf);
   % Nonnegative components summing to the total (A10 contract).
   testCase.verifyEqual(returned_rainf + returned_snowf, ppt);
   testCase.verifyTrue(all(returned_rainf >= 0 & returned_snowf >= 0));
end

function test_initforcings_phase_source_missing_split_stays_missing(testCase)
   % A product without split channels exposes honest all-NaN components in
   % 'source' mode instead of fabricating a split (POLICY A10); the
   % default workspace fixture ships ppt only.
   workspace = testCase.TestData.workspace;
   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016);

   [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, returned_rainf, ...
      returned_snowf] = icemodel.surface.initialize_surface_forcings(opts);

   testCase.verifyTrue(all(isnan(returned_rainf)));
   testCase.verifyTrue(all(isnan(returned_snowf)));
end

function test_resolvePrecipPhase_enforces_validity_contract(testCase)
   % The phase-resolution helper enforces the POLICY A10 contract at every
   % finite sample and refuses unknown options and mismatched axes with
   % named errors.
   ppt = [1e-6; 2e-6];
   tair = [270; 275];

   % Components that do not sum to the total refuse.
   testCase.verifyError(@() icemodel.resolvePrecipPhase( ...
      ppt, tair, [1e-6; 1e-6], [0; 2e-6], 'source'), ...
      'icemodel:resolvePrecipPhase:inconsistentSplit');
   % Negative components refuse even when the sum matches.
   testCase.verifyError(@() icemodel.resolvePrecipPhase( ...
      ppt, tair, [-1e-6; 2e-6], [2e-6; 0], 'source'), ...
      'icemodel:resolvePrecipPhase:inconsistentSplit');
   % A negative partial phase refuses even when its counterpart is missing.
   testCase.verifyError(@() icemodel.resolvePrecipPhase( ...
      ppt, tair, [-1e-6; NaN], [NaN; NaN], 'source'), ...
      'icemodel:resolvePrecipPhase:inconsistentSplit');
   testCase.verifyError(@() icemodel.resolvePrecipPhase( ...
      ppt, tair, [NaN; NaN], [-1e-6; NaN], 'source'), ...
      'icemodel:resolvePrecipPhase:inconsistentSplit');
   % A partial phase cannot exceed its finite total.
   testCase.verifyError(@() icemodel.resolvePrecipPhase( ...
      ppt, tair, [2e-6; NaN], [NaN; NaN], 'source'), ...
      'icemodel:resolvePrecipPhase:inconsistentSplit');
   testCase.verifyError(@() icemodel.resolvePrecipPhase( ...
      ppt, tair, [NaN; NaN], [2e-6; NaN], 'source'), ...
      'icemodel:resolvePrecipPhase:inconsistentSplit');
   % Unknown phase options refuse by name.
   testCase.verifyError(@() icemodel.resolvePrecipPhase( ...
      ppt, tair, ppt, zeros(2, 1), 'bogus'), ...
      'icemodel:resolvePrecipPhase:invalidPhaseSource');
   % Mismatched sample axes refuse before any selection.
   testCase.verifyError(@() icemodel.resolvePrecipPhase( ...
      ppt, [270; 275; 280], ppt, zeros(2, 1), 'source'), ...
      'icemodel:resolvePrecipPhase:sizeMismatch');

   % A consistent split passes through bit-identically with missing
   % samples preserved (honest missingness in 'source' mode).
   rainf = [1e-6; NaN];
   snowf = [0; NaN];
   [returned_rainf, returned_snowf] = icemodel.resolvePrecipPhase( ...
      ppt, tair, rainf, snowf, 'source');
   testCase.verifyEqual(returned_rainf, rainf);
   testCase.verifyEqual(returned_snowf, snowf);
   [returned_rainf, returned_snowf] = icemodel.resolvePrecipPhase( ...
      ppt, tair, [0.5e-6; NaN], [NaN; 1e-6], 'source');
   testCase.verifyEqual(returned_rainf, [0.5e-6; NaN]);
   testCase.verifyEqual(returned_snowf, [NaN; 1e-6]);
end

function test_loadmet_standardizes_optional_snow_depth_alias(testCase)
   % LOADMET should expose one canonical snow_depth series when forcing
   % files carry a supported alias such as snowd.

   workspace = testCase.TestData.workspace;
   snow_depth = linspace(0, 0.12, workspace.nsteps)';
   icemodel.test.fixtures.makeSyntheticMetFile(2016, ...
      'sitename', workspace.sitename, ...
      'forcings', workspace.forcings, ...
      'nsteps', workspace.nsteps, ...
      'dt_seconds', workspace.dt_seconds, ...
      'snow_depth', snow_depth, ...
      'metdir', workspace.metdir);

   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016);
   met = icemodel.loadmet(opts);

   testCase.verifyTrue(isvariable('snow_depth', met));
   testCase.verifyEqual(met.snow_depth, snow_depth, 'AbsTol', 1e-12);
end

function test_loadmet_derives_total_precip_from_components(testCase)
   % LOADMET should expose one canonical total-precip series at runtime
   % when a forcing file carries the rainf/snowf split without ppt (or
   % with a placeholder), so sources like MAR never require a restage to
   % become precipitation-ready downstream.

   workspace = testCase.TestData.workspace;
   rainf = linspace(0, 2e-7, workspace.nsteps)';
   snowf = 5e-8 + zeros(workspace.nsteps, 1);
   % One missing component sample must stay missing in the derived total:
   % the derivation never fabricates a total from half the split.
   rainf(3) = NaN;
   icemodel.test.fixtures.makeSyntheticMetFile(2016, ...
      'sitename', workspace.sitename, ...
      'forcings', workspace.forcings, ...
      'nsteps', workspace.nsteps, ...
      'dt_seconds', workspace.dt_seconds, ...
      'ppt', nan(workspace.nsteps, 1), ...
      'rainf', rainf, ...
      'snowf', snowf, ...
      'metdir', workspace.metdir);

   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016);
   met = icemodel.loadmet(opts);

   testCase.verifyTrue(isvariable('ppt', met));
   expected = rainf + snowf;
   returned = met.ppt;
   testCase.verifyEqual(returned(~isnan(expected)), ...
      expected(~isnan(expected)), 'AbsTol', 1e-15);
   testCase.verifyTrue(isnan(returned(3)));

   % Restore the default fixture: the split channels exist only in this
   % test's file, and a later multi-year load would fail to concatenate
   % timetables with mismatched variables.
   icemodel.test.fixtures.makeSyntheticMetFile(2016, ...
      'sitename', workspace.sitename, ...
      'forcings', workspace.forcings, ...
      'nsteps', workspace.nsteps, ...
      'dt_seconds', workspace.dt_seconds, ...
      'metdir', workspace.metdir);
end

function test_loadmet_exposes_ppt_for_placeholder_split(testCase)
   % A file whose split channels are entirely nonfinite placeholders must
   % still expose the canonical ppt column (all-NaN): downstream
   % contracts require one total-precip channel whenever the split is
   % present, even when nothing was derivable.

   workspace = testCase.TestData.workspace;
   icemodel.test.fixtures.makeSyntheticMetFile(2015, ...
      'sitename', workspace.sitename, ...
      'forcings', workspace.forcings, ...
      'nsteps', workspace.nsteps, ...
      'dt_seconds', workspace.dt_seconds, ...
      'ppt', nan(workspace.nsteps, 1), ...
      'rainf', nan(workspace.nsteps, 1), ...
      'snowf', nan(workspace.nsteps, 1), ...
      'metdir', workspace.metdir);

   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2015);
   met = icemodel.loadmet(opts);

   testCase.verifyTrue(isvariable('ppt', met));
   testCase.verifyTrue(all(isnan(met.ppt)));

   % Restore the default fixture for the same concatenation reason as the
   % derived-total test above.
   icemodel.test.fixtures.makeSyntheticMetFile(2015, ...
      'sitename', workspace.sitename, ...
      'forcings', workspace.forcings, ...
      'nsteps', workspace.nsteps, ...
      'dt_seconds', workspace.dt_seconds, ...
      'metdir', workspace.metdir);
end

function test_initforcings_returns_optional_snow_depth_vector(testCase)
   % icemodel.surface.initialize_surface_forcings should expose the canonical
   % forcing snow-depth vector so the model mains can choose whether to use it.

   workspace = testCase.TestData.workspace;
   snow_depth = 0.05 + zeros(workspace.nsteps, 1);
   icemodel.test.fixtures.makeSyntheticMetFile(2016, ...
      'sitename', workspace.sitename, ...
      'forcings', workspace.forcings, ...
      'nsteps', workspace.nsteps, ...
      'dt_seconds', workspace.dt_seconds, ...
      'snow_depth', snow_depth, ...
      'metdir', workspace.metdir);

   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016);

   [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, met_snow_depth] = ...
      icemodel.surface.initialize_surface_forcings(opts);

   testCase.verifyEqual(met_snow_depth, snow_depth, 'AbsTol', 1e-12);
end

function test_loadmet_swaps_inline_modis_from_metfile(testCase)
   % Inline MODIS-style data embedded in the met file should override the
   % requested output variable when userdata is selected.

   workspace = testCase.TestData.workspace;
   [met_src, ~] = icemodel.test.fixtures.makeSyntheticMetFile(2016, ...
      'sitename', workspace.sitename, ...
      'forcings', workspace.forcings, ...
      'nsteps', workspace.nsteps, ...
      'dt_seconds', workspace.dt_seconds, ...
      'include_modis', true, ...
      'metdir', workspace.metdir);

   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016, userdata='modis', uservars='albedo');
   met = icemodel.loadmet(opts);

   testCase.verifyEqual(met.albedo, met_src.MODIS, 'AbsTol', 1e-12);
end

function test_loadmet_swaps_external_userdata_file(testCase)
   % External userdata files should override the selected met variable when
   % the requested timetable contains the expected replacement field.

   workspace = testCase.TestData.workspace;
   opts_base = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016);
   met_base = icemodel.loadmet(opts_base);

   userdata_values = min(max(met_base.albedo - 0.08, 0.05), 0.95);
   icemodel.test.fixtures.writeSyntheticUserdataFile( ...
      workspace.userdatadir, 2016, ...
      'sitename', workspace.sitename, ...
      'userdata', 'modis', ...
      'varname', 'modis', ...
      'Time', met_base.Time, ...
      'values', userdata_values);

   opts_swap = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016, userdata='modis', uservars='albedo');
   met_swap = icemodel.loadmet(opts_swap);

   testCase.verifyEqual(met_swap.albedo, userdata_values, 'AbsTol', 1e-12);
end

function test_loadmet_swaps_external_userdata_window_file(testCase)
   % A full-period WINDOW userdata file (writeuserdata naming="window":
   % <site>_<source>_<YYYYMMDD>_<YYYYMMDD>.mat) is resolved and used for a run
   % year its encoded period brackets - the loader prefers it over a per-year
   % file and retimes the full-period Data onto the run-year met axis.

   workspace = testCase.TestData.workspace;
   opts_base = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016);
   met_base = icemodel.loadmet(opts_base);

   userdata_values = min(max(met_base.albedo - 0.08, 0.05), 0.95);
   % Window Data spans 2015-2017 and contains the 2016 met timestamps so the
   % linear retime recovers the values exactly at the run year.
   t = [met_base.Time(1) - calmonths(6); met_base.Time; ...
      met_base.Time(end) + calmonths(6)];
   v = [userdata_values(1); userdata_values; userdata_values(end)];
   Data = timetable(v, 'RowTimes', t, 'VariableNames', {'modis'});
   fname = sprintf('%s_modis_%s_%s.mat', workspace.sitename, ...
      char(min(t), 'yyyyMMdd'), char(max(t), 'yyyyMMdd'));
   save(fullfile(workspace.userdatadir, fname), 'Data');

   opts_swap = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016, userdata='modis', uservars='albedo');
   met_swap = icemodel.loadmet(opts_swap);

   testCase.verifyEqual(met_swap.albedo, userdata_values, 'AbsTol', 1e-12);
end

function test_loadmet_explicit_userdata_path_selects_native_variant(testCase)
   % A manifest-selected native Data path must win over the cadence-blind legacy
   % hourly filename when both share one site/source/window identity.
   workspace = testCase.TestData.workspace;
   opts_base = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016);
   met_base = icemodel.loadmet(opts_base);

   hourly_file = fullfile(workspace.userdatadir, ...
      sprintf('%s_modis_20160101_20160101.mat', workspace.sitename));
   Data = timetable(repmat(0.7, height(met_base), 1), ...
      'RowTimes', met_base.Time, 'VariableNames', {'modis'});
   save(hourly_file, 'Data');

   native_time = (met_base.Time(1):minutes(30):met_base.Time(end))';
   Data = timetable(repmat(0.3, numel(native_time), 1), ...
      'RowTimes', native_time, 'VariableNames', {'modis'});
   native_file = fullfile(workspace.userdatadir, ...
      sprintf('%s_modis_20160101_20160101_30m.mat', workspace.sitename));
   save(native_file, 'Data');

   opts_swap = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016, userdata='modis', uservars='albedo');
   missing_year = fullfile(workspace.userdatadir, ...
      sprintf('%s_modis_2015_30m.mat', workspace.sitename));
   opts_swap = icemodel.resetopts(opts_swap, ...
      'userdatafname', {missing_year, native_file});
   met_swap = icemodel.loadmet(opts_swap);

   testCase.verifyEqual(met_swap.albedo, ...
      repmat(0.3, height(met_swap), 1), 'AbsTol', 1e-12);
end

function test_loadmet_explicit_userdata_requires_covering_support(testCase)
   % An existing explicit artifact outside the requested support must fail
   % instead of retiming/extrapolating the wrong Data variant.
   workspace = testCase.TestData.workspace;
   Time = (datetime(2015, 1, 1, 'TimeZone', 'UTC'):hours(1): ...
      datetime(2015, 1, 1, 23, 0, 0, 'TimeZone', 'UTC'))';
   Data = timetable(Time, repmat(0.4, numel(Time), 1), ...
      'VariableNames', {'modis'});
   outside_file = fullfile(workspace.userdatadir, ...
      sprintf('%s_modis_2015_30m.mat', workspace.sitename));
   save(outside_file, 'Data');

   opts_swap = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016, userdata='modis', uservars='albedo');
   opts_swap = icemodel.resetopts(opts_swap, ...
      'userdatafname', {outside_file});

   testCase.verifyError(@() icemodel.loadmet(opts_swap), ...
      'icemodel:loadmet:explicitUserdataCoverage');
end

function test_loadmet_resolves_met_file_in_source_subfolder(testCase)
   % The runtime resolves a met file staged in the per-source subfolder
   % met/<forcings>/ (the writemet staging layout). Staging ONLY in the
   % subfolder (no flat file for this year) proves subfolder-first resolution.

   workspace = testCase.TestData.workspace;
   metsub = fullfile(workspace.metdir, workspace.forcings);
   if ~isfolder(metsub)
      mkdir(metsub);
   end
   icemodel.test.fixtures.makeSyntheticMetFile(2017, ...
      'sitename', workspace.sitename, ...
      'forcings', workspace.forcings, ...
      'nsteps', workspace.nsteps, ...
      'dt_seconds', workspace.dt_seconds, ...
      'metdir', metsub);

   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2017);
   met = icemodel.loadmet(opts);

   testCase.verifyEqual(unique(year(met.Time))', 2017);
   testCase.verifyEqual(height(met), workspace.nsteps);
end

function test_loadmet_prefers_subfolder_met_over_flat(testCase)
   % When the same-named met file exists both flat and in the per-source
   % subfolder, the subfolder copy wins (subfolder-first is the single rule).
   % The setup() already wrote a flat 2016 file; stage a subfolder 2016 file
   % whose tair is offset so the resolved data is unambiguously the subfolder.

   workspace = testCase.TestData.workspace;
   metsub = fullfile(workspace.metdir, workspace.forcings);
   if ~isfolder(metsub)
      mkdir(metsub);
   end
   [met_sub, ~] = icemodel.test.fixtures.makeSyntheticMetFile(2016, ...
      'sitename', workspace.sitename, ...
      'forcings', workspace.forcings, ...
      'nsteps', workspace.nsteps, ...
      'dt_seconds', workspace.dt_seconds, ...
      'metdir', metsub);
   met_sub.tair = met_sub.tair + 5;   % distinguish from the flat fixture
   S.met = met_sub;
   metname = sprintf('met_%s_%s_2016_1hr.mat', ...
      workspace.sitename, workspace.forcings);
   save(fullfile(metsub, metname), '-struct', 'S');

   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016);
   met = icemodel.loadmet(opts);

   testCase.verifyEqual(met.tair, met_sub.tair, 'AbsTol', 1e-12);
end

function test_loadmet_swaps_userdata_in_source_subfolder(testCase)
   % An external userdata file staged in userdata/<source>/ (the writeuserdata
   % staging layout) is resolved by loadmet, proving subfolder-first userdata
   % resolution alongside the flat path covered above.

   workspace = testCase.TestData.workspace;
   opts_base = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016);
   met_base = icemodel.loadmet(opts_base);

   userdata_values = min(max(met_base.albedo - 0.08, 0.05), 0.95);
   userdatasub = fullfile(workspace.userdatadir, 'modis');
   if ~isfolder(userdatasub)
      mkdir(userdatasub);
   end
   icemodel.test.fixtures.writeSyntheticUserdataFile( ...
      userdatasub, 2016, ...
      'sitename', workspace.sitename, ...
      'userdata', 'modis', ...
      'varname', 'modis', ...
      'Time', met_base.Time, ...
      'values', userdata_values);

   opts_swap = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016, userdata='modis', uservars='albedo');
   met_swap = icemodel.loadmet(opts_swap);

   testCase.verifyEqual(met_swap.albedo, userdata_values, 'AbsTol', 1e-12);
end

function test_loadmet_prefers_external_met_swap_source(testCase)
   % Swap sources should prefer staged met files over legacy userdata files so
   % RCM met artifacts can fill missing native channels without duplicate Data
   % files. The target run is 15 min while the source met is hourly, matching
   % the PROMICE-plus-RCM development workflow.

   workspace = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
      configure=true, nsteps=8, dt_seconds=900);
   cleanup = onCleanup(@() ...
      icemodel.test.fixtures.cleanupSyntheticWorkspace(workspace));

   opts_base = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016);
   met_base = icemodel.loadmet(opts_base);

   metsub = fullfile(workspace.metdir, 'modis');
   if ~isfolder(metsub)
      mkdir(metsub);
   end
   [met_source, source_file] = icemodel.test.fixtures.makeSyntheticMetFile( ...
      2016, sitename=workspace.sitename, forcings='modis', nsteps=3, ...
      dt_seconds=3600, metdir=metsub);
   met_source.albedo = 0.35 + (0:height(met_source)-1)' * 0.01;
   met = met_source;
   save(source_file, 'met');

   userdata_values = 0.91 + zeros(height(met_base), 1);
   icemodel.test.fixtures.writeSyntheticUserdataFile( ...
      workspace.userdatadir, 2016, ...
      'sitename', workspace.sitename, ...
      'userdata', 'modis', ...
      'varname', 'modis', ...
      'Time', met_base.Time, ...
      'values', userdata_values);

   opts_swap = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016, userdata='modis', uservars='albedo');
   met_swap = icemodel.loadmet(opts_swap);

   expected = retime(met_source, met_base.Time, 'linear');
   testCase.verifyEqual(met_swap.albedo, expected.albedo, 'AbsTol', 1e-12);
   testCase.verifyNotEqual(met_swap.albedo, userdata_values);

   clear cleanup
end

function test_loadmet_supports_native_30m_swap_source_and_fallbacks(testCase)
   % A native 30-minute run must prefer exact 30m swap met, then retain the
   % established compatible 15-minute and hourly source fallbacks.
   workspace = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
      configure=true, nsteps=8, dt_seconds=1800);
   cleanup = onCleanup(@() ...
      icemodel.test.fixtures.cleanupSyntheticWorkspace(workspace));
   metsub = fullfile(workspace.metdir, 'modis');
   mkdir(metsub);

   % Distinct constant albedos make the selected cadence observable after the
   % source timetable is retimed onto the target 30-minute run.
   [met_30m, file_30m] = icemodel.test.fixtures.makeSyntheticMetFile( ...
      2016, sitename=workspace.sitename, forcings='modis', nsteps=8, ...
      dt_seconds=1800, metdir=metsub);
   met_30m.albedo(:) = 0.30;
   met = met_30m;
   save(file_30m, 'met');
   [met_15m, file_15m] = icemodel.test.fixtures.makeSyntheticMetFile( ...
      2016, sitename=workspace.sitename, forcings='modis', nsteps=15, ...
      dt_seconds=900, metdir=metsub);
   met_15m.albedo(:) = 0.15;
   met = met_15m;
   save(file_15m, 'met');
   [met_hourly, file_hourly] = icemodel.test.fixtures.makeSyntheticMetFile( ...
      2016, sitename=workspace.sitename, forcings='modis', nsteps=5, ...
      dt_seconds=3600, metdir=metsub);
   met_hourly.albedo(:) = 0.60;
   met = met_hourly;
   save(file_hourly, 'met');

   opts_swap = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016, userdata='modis', uservars='albedo');
   selected = icemodel.loadmet(opts_swap);
   testCase.verifyEqual(selected.albedo, ...
      repmat(0.30, height(selected), 1), 'AbsTol', 1e-12);

   % Removing exact support exposes the documented finer then hourly fallback
   % order without changing the target run cadence.
   delete(file_30m);
   selected = icemodel.loadmet(opts_swap);
   testCase.verifyEqual(selected.albedo, ...
      repmat(0.15, height(selected), 1), 'AbsTol', 1e-12);
   delete(file_15m);
   selected = icemodel.loadmet(opts_swap);
   testCase.verifyEqual(selected.albedo, ...
      repmat(0.60, height(selected), 1), 'AbsTol', 1e-12);

   clear cleanup
end

function test_loadmet_errors_when_only_daily_userdata_swap_available(testCase)
   % Daily observations are not safe met-channel swap sources. If no matching
   % met file exists, loadmet should fail clearly instead of interpolating a
   % daily userdata fallback into subdaily forcing.

   workspace = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
      configure=true, nsteps=24, dt_seconds=3600);
   cleanup = onCleanup(@() ...
      icemodel.test.fixtures.cleanupSyntheticWorkspace(workspace));

   Time = [datetime(2016, 1, 1, 'TimeZone', 'UTC'); ...
      datetime(2016, 1, 2, 'TimeZone', 'UTC')];
   modis = [0.45; 0.46];
   Data = timetable(Time, modis);
   save(fullfile(workspace.userdatadir, 'kanm_modis_2016.mat'), 'Data');

   opts_swap = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016, userdata='modis', uservars='albedo');

   testCase.verifyError(@() icemodel.loadmet(opts_swap), ...
      'icemodel:loadmet:dailySwapData');

   clear cleanup
end

function test_loadmet_errors_when_userdata_file_lacks_Data(testCase)
   % Corrupt userdata files must raise an error instead of falling back to
   % the met data.

   workspace = testCase.TestData.workspace;
   filepath = fullfile(workspace.userdatadir, 'kanm_modis_2016.mat');
   warn_state = warning('off', 'MATLAB:load:variableNotFound');
   cleanup = onCleanup(@() warning(warn_state));
   bogus = 1;
   save(filepath, 'bogus');

   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016, userdata='modis', uservars='albedo');

   try
      icemodel.loadmet(opts);
      testCase.verifyFail('expected loadmet to fail when userdata lacks Data');
   catch ME
      testCase.verifyTrue(contains(ME.message, ...
         'userdata file does not contain timetable "Data"'));
   end
   clear cleanup
end

function test_processmet_supports_native_and_hourly_cadence(testCase)
   % PROCESSMET should preserve native cadence when requested and collapse
   % to hourly outputs when the hourly contract is requested instead.

   [met_native, ~] = icemodel.test.fixtures.makeSyntheticMetFile(2016, ...
      'nsteps', 8, 'dt_seconds', 900);

   native = icemodel.processmet(met_native, newTimeStep="native");
   hourly = icemodel.processmet(met_native, newTimeStep="hourly");

   testCase.verifyEqual(height(native), 8);
   testCase.verifyEqual(height(hourly), 2);
   testCase.verifyTrue(all(ismember( ...
      {'tsfc', 'swu', 'swn', 'lwu', 'lwn', 'netr'}, ...
      string(hourly.Properties.VariableNames))));
end

function test_loadresults_defaults_to_output_years(testCase)
   % LOADRESULTS should prefer opts.output_years when the saved run spans
   % spinup and retained output years.

   workspace = testCase.TestData.workspace;
   pathoutput = fullfile(workspace.outputdir, 'kanm', 'skinmodel', ...
      'loadresults');
   if exist(pathoutput, 'dir') ~= 7
      mkdir(pathoutput);
   end
   if exist(fullfile(pathoutput, '2016'), 'dir') ~= 7
      mkdir(fullfile(pathoutput, '2016'));
   end
   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', [2015 2016], ...
      testname='loadresults', saveflag=true, n_spinup_years=1, ...
      output_profile='minimal', solver=1, dt=3600, pathoutput=pathoutput);
   [~, ~, opts] = icemodel.test.helpers.runSmbModel(opts);

   [ice1, ~, met] = icemodel.loadresults(opts);

   testCase.verifyEqual(unique(year(ice1.Time))', 2016);
   testCase.verifyEqual(unique(year(met.Time))', 2016);
end

function test_loadresults_preserves_diagnostic_output_profile_fields(testCase)
   % Saved diagnostic-profile outputs should reload with the accepted
   % solver/THF debugging scalars intact.

   workspace = testCase.TestData.workspace;
   pathoutput = fullfile(workspace.outputdir, 'kanm', 'icemodel', ...
      'loadresults_diagnostic');
   if exist(pathoutput, 'dir') ~= 7
      mkdir(pathoutput);
   end
   if exist(fullfile(pathoutput, '2016'), 'dir') ~= 7
      mkdir(fullfile(pathoutput, '2016'));
   end

   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'icemodel', 2016, saveflag=true, ...
      output_profile='diagnostic', solver=1, seb_solver=2, ...
      turbulent_flux_scheme='monin_obukhov', z0_ice=0.02, ...
      pathoutput=pathoutput, testname='loadresults_diagnostic');
   [~, ~, opts] = icemodel.test.helpers.runSmbModel(opts);

   [ice1, ~, met] = icemodel.loadresults(opts);

   testCase.verifyTrue(all(ismember( ...
      {'n_subfail', 'ea_atm', 'br_coefs_gamma', 'br_coefs_b1_num', ...
      'br_coefs_b2_num', 'hv_atm', 'ro_sfc', ...
      'thf_es_sfc', 'thf_z0m', 'thf_z0h', 'thf_z0q', 'thf_u_star', ...
      'thf_L', 'thf_Re', 'thf_numiter'}, ice1.Properties.VariableNames)));
   testCase.verifyTrue(all(isfinite(ice1.ro_sfc)));
   testCase.verifyTrue(all(isfinite(ice1.thf_u_star)));
   testCase.verifyTrue(all(isfinite(ice1.thf_L)));
   testCase.verifyTrue(all(isfinite(ice1.thf_numiter)));

   % Diagnostic saves must keep every step-ledger channel and its
   % signed storage, remesh, geometry, and phase-aware vapor identities.
   budget_fields = icemodel.namelists.budgetoutputs();
   testCase.verifyTrue(all(ismember( ...
      budget_fields, ice1.Properties.VariableNames)));
   testCase.verifyTrue(all(isfinite(ice1{:, budget_fields}), 'all'));

   % Every absolute-gross channel is nonnegative and bounds its signed
   % forcing-step net after accepted substeps and remesh events are accumulated.
   gross_fields = budget_fields(contains(budget_fields, '_gross_'));
   testCase.verifyGreaterThanOrEqual(ice1{:, gross_fields}, 0);

   solid_delta = ice1.mass_budget_solid_end_mwe ...
      - ice1.mass_budget_solid_start_mwe;
   liquid_delta = ice1.mass_budget_liquid_end_mwe ...
      - ice1.mass_budget_liquid_start_mwe;
   testCase.verifyEqual(solid_delta, ...
      ice1.mass_budget_phase_solid_mwe ...
      + ice1.mass_budget_vapor_solid_mwe ...
      + ice1.mass_budget_remesh_solid_mwe, 'AbsTol', 1e-10);
   testCase.verifyEqual(liquid_delta, ...
      ice1.mass_budget_phase_liquid_mwe ...
      + ice1.mass_budget_vapor_liquid_mwe ...
      + ice1.mass_budget_remesh_liquid_mwe, 'AbsTol', 1e-10);
   testCase.verifyGreaterThanOrEqual( ...
      ice1.mass_budget_solid_storage_gross_mwe + 1e-12, ...
      abs(solid_delta));
   testCase.verifyGreaterThanOrEqual( ...
      ice1.mass_budget_liquid_storage_gross_mwe + 1e-12, ...
      abs(liquid_delta));

   signed_gross_pairs = { ...
      'mass_budget_phase_solid_mwe', ...
      'mass_budget_phase_solid_gross_mwe'; ...
      'mass_budget_phase_liquid_mwe', ...
      'mass_budget_phase_liquid_gross_mwe'; ...
      'mass_budget_vapor_solid_mwe', ...
      'mass_budget_vapor_solid_gross_mwe'; ...
      'mass_budget_vapor_liquid_mwe', ...
      'mass_budget_vapor_liquid_gross_mwe'; ...
      'mass_budget_remesh_solid_mwe', ...
      'mass_budget_remesh_solid_gross_mwe'; ...
      'mass_budget_remesh_liquid_mwe', ...
      'mass_budget_remesh_liquid_gross_mwe'; ...
      'mass_budget_vapor_potential_j_m2', ...
      'mass_budget_vapor_potential_gross_j_m2'; ...
      'mass_budget_unapplied_vapor_j_m2', ...
      'mass_budget_unapplied_vapor_gross_j_m2'};
   for n = 1:size(signed_gross_pairs, 1)
      % A single tolerance only absorbs binary roundoff in the comparison.
      signed_name = signed_gross_pairs{n, 1};
      gross_name = signed_gross_pairs{n, 2};
      testCase.verifyGreaterThanOrEqual( ...
         ice1.(gross_name) + 1e-8, abs(ice1.(signed_name)));
   end
   testCase.verifyEqual(ice1.mass_budget_remesh_solid_mwe, ...
      ice1.mass_budget_cloned_bottom_solid_mwe ...
      - ice1.mass_budget_merge_export_solid_mwe, 'AbsTol', 1e-10);
   testCase.verifyEqual(ice1.mass_budget_remesh_liquid_mwe, ...
      ice1.mass_budget_cloned_bottom_liquid_mwe ...
      - ice1.mass_budget_merge_export_liquid_mwe, 'AbsTol', 1e-10);
   testCase.verifyEqual(ice1.mass_budget_remesh_enthalpy_j_m2, ...
      ice1.mass_budget_cloned_bottom_enthalpy_j_m2 ...
      - ice1.mass_budget_merge_export_enthalpy_j_m2, 'AbsTol', 1e-6);

   [Ls, Lv, ro_liq] = icemodel.physicalConstant('Ls', 'Lv', 'ro_liq');
   vapor_accounted = ro_liq * ( ...
      Ls * ice1.mass_budget_vapor_solid_mwe ...
      + Lv * ice1.mass_budget_vapor_liquid_mwe ...
      + Lv * ice1.mass_budget_condensation_overflow_mwe) ...
      + ice1.mass_budget_unapplied_vapor_j_m2;
   testCase.verifyEqual(ice1.mass_budget_vapor_potential_j_m2, ...
      vapor_accounted, 'AbsTol', 2e-5);
   testCase.verifyEqual(ice1.mass_budget_top_deletion_height_m, ...
      opts.dz_thermal * ice1.mass_budget_top_deletion_count, ...
      'AbsTol', 1e-12);
   testCase.verifyEqual(unique(year(met.Time))', 2016);
end

function test_postprocess_explicit_met_return_is_hourly(testCase)
   % The explicit met output from POSTPROCESS should stay aligned with the
   % hourly postprocessed model outputs.

   localws = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
      configure=true, nsteps=96, dt_seconds=900);
   cleanup = onCleanup(@() icemodel.test.fixtures.cleanupSyntheticWorkspace( ...
      localws));

   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      localws, 'skinmodel', 2016, dt=900, solver=1, testname='postprocess');
   [ice1_raw, ice2_raw, opts] = icemodel.test.helpers.runSmbModel(opts);
   met = icemodel.loadmet(opts);

   [ice1_pp, ~, met_pp] = icemodel.postprocess(ice1_raw, ice2_raw, opts, ...
      met.swd, met.lwd, met.albedo, met.Time);

   testCase.verifyEqual(height(ice1_pp), 24);
   testCase.verifyEqual(height(met_pp), 24);
   clear cleanup
end

function test_quarter_hour_diagnostic_postprocess_preserves_budget_contract(testCase)
   % Hourly diagnostic output must retain double ledgers and closure identities.

   localws = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
      configure=true, nsteps=8, dt_seconds=900);
   cleanup = onCleanup(@() icemodel.test.fixtures.cleanupSyntheticWorkspace( ...
      localws));
   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      localws, 'icemodel', 2016, dt=900, solver=1, ...
      output_profile='diagnostic', testname='diagnostic_retime');
   [ice1_raw, ice2_raw, opts] = icemodel.test.helpers.runSmbModel(opts);
   met = icemodel.loadmet(opts);

   [ice1, ~] = icemodel.postprocess(ice1_raw, ice2_raw, opts, ...
      met.swd, met.lwd, met.albedo, met.Time);

   % Every ledger stays double even though convergence flags remain single.
   budget_fields = icemodel.namelists.budgetoutputs();
   for n = 1:numel(budget_fields)
      testCase.verifyClass(ice1.(budget_fields{n}), 'double');
   end
   testCase.verifyClass(ice1.Tice_converged, 'single');

   % Sum/first/last hourly aggregation preserves the signed checkpoint closure.
   solid_delta = ice1.mass_budget_solid_end_mwe ...
      - ice1.mass_budget_solid_start_mwe;
   liquid_delta = ice1.mass_budget_liquid_end_mwe ...
      - ice1.mass_budget_liquid_start_mwe;
   testCase.verifyEqual(solid_delta, ...
      ice1.mass_budget_phase_solid_mwe ...
      + ice1.mass_budget_vapor_solid_mwe ...
      + ice1.mass_budget_remesh_solid_mwe, 'AbsTol', 1e-10);
   testCase.verifyEqual(liquid_delta, ...
      ice1.mass_budget_phase_liquid_mwe ...
      + ice1.mass_budget_vapor_liquid_mwe ...
      + ice1.mass_budget_remesh_liquid_mwe, 'AbsTol', 1e-10);
   testCase.verifyGreaterThanOrEqual( ...
      ice1.mass_budget_solid_storage_gross_mwe + 1e-12, ...
      abs(solid_delta));
   testCase.verifyGreaterThanOrEqual( ...
      ice1.mass_budget_liquid_storage_gross_mwe + 1e-12, ...
      abs(liquid_delta));
   clear cleanup
end

function test_postprocess_three_row_partial_bin_aligns_subsurface_and_logical(testCase)
   % A partial native hour must use all available samples for every output.

   localws = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
      configure=true, nsteps=3, dt_seconds=900);
   cleanup = onCleanup(@() icemodel.test.fixtures.cleanupSyntheticWorkspace( ...
      localws));
   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      localws, 'icemodel', 2016, dt=900, solver=1, ...
      output_profile='minimal', testname='partial_hour_retime');
   [ice1_raw, ice2_raw, opts] = icemodel.test.helpers.runSmbModel(opts);
   met = icemodel.loadmet(opts);

   % Seed every layer with values whose three-sample mean and logical union
   % are unambiguous after hourly postprocessing.
   n_layers = size(ice2_raw.Tice, 1);
   ice2_raw.bin_sentinel = repmat([1, 2, 3], n_layers, 1);
   ice2_raw.bin_flag = repmat(logical([0, 1, 0]), n_layers, 1);
   [ice1, ice2] = icemodel.postprocess(ice1_raw, ice2_raw, opts, ...
      met.swd, met.lwd, met.albedo, met.Time);

   % The subsurface mean and logical flag must occupy the same sole surface
   % bin rather than being dropped because a fourth sample is unavailable.
   expected_time = dateshift(met.Time(1), 'start', 'hour');
   testCase.verifyEqual(ice1.Properties.RowTimes, expected_time);
   testCase.verifyEqual(ice2.bin_sentinel, 2 * ones(n_layers, 1));
   testCase.verifyClass(ice1.bin_flag, 'logical');
   testCase.verifyTrue(ice1.bin_flag);
   clear cleanup
end

function test_postprocess_unaligned_bins_align_subsurface_and_logical(testCase)
   % Unaligned samples must follow native hourly boundaries for every output.

   localws = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
      configure=true, nsteps=4, dt_seconds=900);
   cleanup = onCleanup(@() icemodel.test.fixtures.cleanupSyntheticWorkspace( ...
      localws));
   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      localws, 'icemodel', 2016, dt=900, solver=1, ...
      output_profile='minimal', testname='unaligned_hour_retime');
   [ice1_raw, ice2_raw, opts] = icemodel.test.helpers.runSmbModel(opts);
   met = icemodel.loadmet(opts);

   % Shift the first three samples into one partial hour and the final sample
   % two hours later. The native result includes one intervening empty bin.
   shifted_time = met.Time(1) + minutes([15; 30; 45; 120]);
   n_layers = size(ice2_raw.Tice, 1);
   ice2_raw.bin_sentinel = repmat([1, 2, 3, 100], n_layers, 1);
   ice2_raw.bin_flag = repmat(logical([0, 1, 0, 1]), n_layers, 1);
   [ice1, ice2] = icemodel.postprocess(ice1_raw, ice2_raw, opts, ...
      met.swd, met.lwd, met.albedo, shifted_time);

   % Surface labels, subsurface means, and logical unions must all preserve the
   % first partial bin, the empty hour, and the final one-sample bin.
   expected_time = dateshift(shifted_time(1), 'start', 'hour') ...
      + hours((0:2)');
   testCase.verifyEqual(ice1.Properties.RowTimes, expected_time);
   testCase.verifyEqual(ice2.bin_sentinel(:, [1, 3]), ...
      repmat([2, 100], n_layers, 1));
   testCase.verifyTrue(all(isnan(ice2.bin_sentinel(:, 2)), 'all'));
   testCase.verifyClass(ice1.bin_flag, 'logical');
   testCase.verifyEqual(ice1.bin_flag, [true; false; true]);
   clear cleanup
end

function test_retimeHourlyFixedStep_matches_legacy_hourly_mean(testCase)
   % The fixed-step hourly aggregation helper should reproduce the legacy
   % timetable RETIME path exactly on aligned 15-minute model output.

   localws = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
      configure=true, nsteps=96, dt_seconds=900);
   cleanup = onCleanup(@() icemodel.test.fixtures.cleanupSyntheticWorkspace( ...
      localws));

   % Run one synthetic quarter-hour case and convert the raw surface output
   % into the same timetable form POSTPROCESS retimes.
   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      localws, 'skinmodel', 2016, dt=900, solver=1, ...
      testname='postprocess_retime_equivalence');
   [ice1_raw, ~, opts] = icemodel.test.helpers.runSmbModel(opts);
   met = icemodel.loadmet(opts);
   ice1_tt = rawIce1ToTimetable(ice1_raw, met.Time);

   % Compare the legacy timetable RETIME result to the fixed-step helper.
   legacy = legacyHourlyMean(ice1_tt);
   fixed = icemodel.retimeHourlyFixedStep(ice1_tt);

   testCase.verifyEqual(fixed.Time, legacy.Time);
   variables = legacy.Properties.VariableNames;
   for n = 1:numel(variables)
      % Compare each channel independently so table concatenation cannot hide
      % a double-to-single conversion in either result.
      name = variables{n};
      testCase.verifyClass(fixed.(name), class(legacy.(name)));
      testCase.verifyEqual(fixed.(name), legacy.(name), 'AbsTol', 5e-5);
   end
   clear cleanup
end

function test_retimeHourlyFixedStep_uses_cumulative_endpoint_state(testCase)
   % Cumulative diagnostics retain the final native state; ordinary channels
   % retain the established hourly mean.

   time = datetime(2016, 1, 1, 0, 0, 0, 'TimeZone', 'UTC') ...
      + minutes((0:7)' * 15);
   cumulative = [1:4, 11:14]';
   ordinary = [2:2:8, 20:2:26]';
   input = timetable(cumulative, cumulative + 100, cumulative + 200, ...
      cumulative + 300, ordinary, 'RowTimes', time, ...
      'VariableNames', [icemodel.namelists.cumulativeoutputs(), ...
      {'ordinary'}]);

   [hourly, bin_start, bin_end] = ...
      icemodel.retimeHourlyFixedStep(input);

   % The aligned 15-minute grid exercises the fixed-array branch and exposes
   % the same four-sample bounds consumed by subsurface postprocessing.
   testCase.verifyEqual(bin_start, [1; 5]);
   testCase.verifyEqual(bin_end, [4; 8]);
   testCase.verifyEqual(icemodel.namelists.cumulativeoutputs(), ...
      {'melt', 'runoff', 'freeze', 'dlayer'});
   testCase.verifyEqual(hourly{:, ...
      icemodel.namelists.cumulativeoutputs()}, ...
      [4, 104, 204, 304; 14, 114, 214, 314]);
   testCase.verifyEqual(hourly.ordinary, [5; 23]);
end

function test_retimeHourlyFixedStep_preserves_irregular_native_bins(testCase)
   % An aligned row count must not trigger fixed aggregation when gaps make the
   % native grid irregular; cumulative channels still retain each bin endpoint.

   time = datetime(2016, 1, 1, 0, 0, 0, 'TimeZone', 'UTC') ...
      + minutes([0; 15; 45; 60; 75; 180; 195; 210]);
   cumulative = transpose(1:8);
   ordinary = 10 * cumulative;
   input = timetable(cumulative, ordinary, 'RowTimes', time, ...
      'VariableNames', {'melt', 'ordinary'});

   [hourly, bin_start, bin_end] = ...
      icemodel.retimeHourlyFixedStep(input);

   % Native hourly labels and raw bounds retain every irregular sample without
   % allowing the fixed four-row reshape to cross either hour boundary.
   expected_time = datetime(2016, 1, 1, (0:3)', 0, 0, ...
      'TimeZone', 'UTC');
   testCase.verifyEqual(hourly.Properties.RowTimes, expected_time);
   testCase.verifyEqual(bin_start, [1; 4; 0; 6]);
   testCase.verifyEqual(bin_end, [3; 5; 0; 8]);
   testCase.verifyEqual(hourly.melt([1; 2; 4]), [3; 5; 8]);
   testCase.verifyTrue(isnan(hourly.melt(3)));
   testCase.verifyEqual(hourly.ordinary([1; 2; 4]), [20; 45; 70]);
   testCase.verifyTrue(isnan(hourly.ordinary(3)));
end

function test_retimeHourlyFixedStep_sums_increment_channels(testCase)
   % Every df_ channel holds one forcing step's change, so aggregation sums it.
   % Averaging one would divide it by the samples per bin, turning an hourly
   % overflow depth into a quarter of the water that ran off.

   % Aligned fixed-array path: four unit increments become one hourly four.
   time = datetime(2016, 1, 1, 0, 0, 0, 'TimeZone', 'UTC') ...
      + minutes((0:7)' * 15);
   df_rof = ones(8, 1);
   ordinary = ones(8, 1);
   hourly = icemodel.retimeHourlyFixedStep( ...
      timetable(df_rof, ordinary, 'RowTimes', time));
   testCase.verifyEqual(hourly.df_rof, [4; 4]);
   testCase.verifyEqual(hourly.ordinary, [1; 1]);

   % Native irregular path must apply the same class, not the default mean.
   time = datetime(2016, 1, 1, 0, 0, 0, 'TimeZone', 'UTC') ...
      + minutes([0; 15; 45; 60]);
   df_rof = ones(4, 1);
   ordinary = ones(4, 1);
   hourly = icemodel.retimeHourlyFixedStep( ...
      timetable(df_rof, ordinary, 'RowTimes', time));
   testCase.verifyEqual(hourly.df_rof, [3; 1]);
   testCase.verifyEqual(hourly.ordinary, [1; 1]);
end

function test_retimeHourlyFixedStep_matches_native_missing_values(testCase)
   % Fixed-array aggregation must retain RETIME's nonmissing mean, sum, and
   % endpoint selection rather than treating a missing edge as the bin value.

   time = datetime(2016, 1, 1, 0, 0, 0, 'TimeZone', 'UTC') ...
      + minutes((0:3)' * 15);
   values = [NaN; 2; 3; NaN];
   input = timetable(values, values, values, values, values, ...
      'RowTimes', time, 'VariableNames', {'ordinary', 'melt', ...
      'mass_budget_solid_start_mwe', 'mass_budget_solid_end_mwe', ...
      'mass_budget_phase_solid_mwe'});

   hourly = icemodel.retimeHourlyFixedStep(input);

   testCase.verifyEqual(hourly.ordinary, 2.5);
   testCase.verifyEqual(hourly.melt, 3);
   testCase.verifyEqual(hourly.mass_budget_solid_start_mwe, 2);
   testCase.verifyEqual(hourly.mass_budget_solid_end_mwe, 3);
   testCase.verifyEqual(hourly.mass_budget_phase_solid_mwe, 5);
end

function ice1_tt = rawIce1ToTimetable(ice1_raw, time)
   %RAWICE1TOTIMETABLE Convert raw model output to POSTPROCESS timetable form.

   % Match POSTPROCESS by converting logical convergence flags to single
   % before building the timetable-backed hourly series.
   if isfield(ice1_raw, 'Tice_converged')
      ice1_raw.Tice_converged = single(ice1_raw.Tice_converged);
   end
   if isfield(ice1_raw, 'Tsfc_converged')
      ice1_raw.Tsfc_converged = single(ice1_raw.Tsfc_converged);
   end

   % Convert the raw struct into the exact timetable shape retimed by
   % POSTPROCESS.
   time.TimeZone = 'UTC';
   ice1_tt = struct2table(ice1_raw);
   ice1_tt = table2timetable(ice1_tt, 'RowTimes', time);
end

function ice1_hourly = legacyHourlyMean(ice1_tt)
   %LEGACYHOURLYMEAN Reproduce the timetable RETIME branch from POSTPROCESS.

   % Match the legacy path exactly, including the leap-day removal.
   ice1_hourly = retime(ice1_tt, 'hourly', 'mean');
   ice1_hourly = ice1_hourly(~(month(ice1_hourly.Time) == 2 ...
      & day(ice1_hourly.Time) == 29), :);
end

function [workspace, opts, met] = makeCodegenPromiceFixture()
   %MAKECODEGENPROMICEFIXTURE Build one canonical pre-codegen forcing input.
   workspace = icemodel.test.fixtures.makeSyntheticWorkspace( ...
      2015, configure=false, nsteps=96, dt_seconds=900);
   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2015);
   opts.forcings = 'promice_filled';
   opts.userdata = 'promice_filled';
   opts.readiness_file = fullfile(workspace.metdir, 'readiness.csv');
   opts.report_inputs_file = fullfile(workspace.metdir, ...
      'report-inputs.json');
   opts.startdate = datetime(2015, 1, 1, 0, 0, 0, 'TimeZone', 'UTC');
   opts.enddate = datetime(2015, 1, 1, 23, 45, 0, 'TimeZone', 'UTC');

   % Keep the bookkeeping ledger explicitly ready; the test proves it cannot
   % substitute for artifact policy and per-channel provenance.
   ledger = table(string(opts.sitename), 2015, "ready", "ready", ...
      'VariableNames', {'site', 'year', ...
      'verdict_icemodel', 'verdict_snowmodel'});
   writetable(ledger, opts.readiness_file)

   % Promote the native synthetic file into the exact filled-product payload
   % shape consumed by the MATLAB verifier before generated execution.
   loaded = load(opts.metfname{1}, 'met');
   met = loaded.met;
   met.Time = opts.startdate + minutes(15) * (0:height(met) - 1)';
   met.swu = met.albedo .* met.swd;
   met.rainf = zeros(height(met), 1);
   met.snowf = zeros(height(met), 1);
   met.boom_height = repmat(2.7, height(met), 1);
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   defaults = icemodel.forcing.reconstruct.setopts();
   channels = unique([defaults.plan_channels, ...
      icemodel.forcing.helpers.precipitationVariables(), ...
      "boom_height"], 'stable');
   for channel = channels
      provenance = repmat(codes.observed, height(met), 1);
      provenance(~isfinite(met.(channel))) = codes.missing;
      met.(channel + "_provenance") = provenance;
   end

   % Stamp the same current identity that the production reconstruction
   % writer records before creating the manifest-pinned file.
   metadata = met.Properties.UserData;
   metadata.site = string(opts.sitename);
   metadata.gapfill_registry = codes;
   metadata.gapfill_seed = 1;
   metadata.gapfill_product = "promice_filled";
   metadata.gapfill_channels = defaults.plan_channels;
   metadata.gapfill_engine_version = string(icemodel.internal.version());
   metadata.gapfill_policy_sha256 = ...
      icemodel.forcing.reconstruct.policySha256();
   metadata.gapfill_donors = string.empty(1, 0);
   met.Properties.UserData = metadata;
   artifact_metadata = icemodel.forcing.helpers.artifactMetadata(met);
   met.Properties.UserData = artifact_metadata;
   filename = fullfile(workspace.metdir, sprintf( ...
      'met_%s_promice_filled_20150101_20150101_15m.mat', opts.sitename));
   save(filename, 'met', 'artifact_metadata')
   opts.metfname = {filename};
   writePromiceRuntimeManifest(opts)
end

function writePromiceRuntimeManifest(opts)
   %WRITEPROMICERUNTIMEMANIFEST Pin synthetic runtime inputs by byte hash.
   paths = [string(opts.readiness_file); string(opts.metfname(:))];
   roles = ["readiness"; repmat("filled", numel(opts.metfname), 1)];
   [data_root, ~] = icemodel.forcing.reconstruct.selectedDataRoot( ...
      string(fileparts(opts.metfname{1})));
   relative_paths = icemodel.verification.setup.relpaths(paths, data_root);
   artifacts = repmat(struct('role', "", 'path', "", 'bytes', 0, ...
      'sha256', ""), numel(paths), 1);
   for k = 1:numel(paths)
      info = dir(paths(k));
      artifacts(k) = struct('role', roles(k), 'path', relative_paths(k), ...
         'bytes', info.bytes, 'sha256', ...
         icemodel.verification.setup.fileSha256(paths(k)));
   end
   manifest = struct('site', string(opts.sitename), ...
      'path_base', "selected_data_root", 'artifacts', artifacts);
   fid = fopen(opts.report_inputs_file, 'w');
   cleaner = onCleanup(@() fclose(fid));
   fprintf(fid, '%s', jsonencode(manifest));
   clear cleaner
end
