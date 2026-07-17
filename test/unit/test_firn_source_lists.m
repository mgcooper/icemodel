function tests = test_firn_source_lists
   %TEST_FIRN_SOURCE_LISTS Verify manifest source-list derivation.
   tests = functiontests(localfunctions);
end

function test_default_staging_roots_ignore_active_data_path(testCase)
   % Firn importers use blank config casenames by default, so the shared root
   % helpers must resolve repo top-level data paths instead of whatever
   % ICEMODEL_DATA_PATH happens to be in the caller's shell.
   old_data = getenv('ICEMODEL_DATA_PATH');
   cleanup = onCleanup(@() setenv('ICEMODEL_DATA_PATH', old_data));
   setenv('ICEMODEL_DATA_PATH', tempname);

   returned_eval = icemodel.verification.helpers.evaluationDataRoot();
   expected_eval = string(fullfile(icemodel.internal.fullpath('data'), ...
      'eval'));
   testCase.verifyEqual(returned_eval, expected_eval);

   returned_input = icemodel.verification.helpers.inputDataRoot();
   expected_input = string(fullfile(icemodel.internal.fullpath('data'), ...
      'input'));
   testCase.verifyEqual(returned_input, expected_input);

   clear cleanup
end

function test_models_are_eval_sources_when_staged(testCase)
   % MAR/MERRA/RACMO can all be comparison targets. MAR/MERRA are also forcing
   % sources only when they carry met files; RACMO remains eval-only.
   colocation = struct();
   colocation.promice = leg(true, "met_kanu_promice.mat");
   colocation.sumup = struct('staged', true);
   colocation.gcnet = leg(true, "met_dye2_gcnet.mat");
   colocation.mar = leg(true, "met_kanu_mar.mat");
   colocation.merra = leg(true, "met_kanu_merra.mat");
   colocation.racmo = struct('staged', true, 'data_files', "kanu_racmo.mat");

   [returned_forcing, returned_eval] = ...
      icemodel.verification.setup.colocationSourceLists(colocation);

   % Native/met-capable sources with met files and readiness stamps should be
   % forcing sources.
   expected_forcing = ["promice"; "gcnet"; "mar3.11"; "merra2"];
   testCase.verifyEqual(returned_forcing, expected_forcing);

   % Staged observation/model legs should be eval sources; RACMO is eval only.
   expected_eval = ["promice_obs"; "sumup_obs"; "gcnet_obs"; ...
      "mar3.11"; "merra2"; "racmo2.3p3"];
   testCase.verifyEqual(returned_eval, expected_eval);
end

function test_evaluation_staging_is_independent_of_runtime_artifacts(testCase)
   % A case-level PROMICE observation bundle remains an evaluation source when
   % a source-selective call intentionally writes no native met/Data files.
   colocation = struct();
   colocation.promice = struct('staged', false, 'eval_staged', true, ...
      'met_files', strings(1, 0), 'data_files', strings(1, 0));

   [returned_forcing, returned_eval] = ...
      icemodel.verification.setup.colocationSourceLists(colocation);

   testCase.verifyEmpty(returned_forcing);
   testCase.verifyEqual(returned_eval, "promice_obs");
end

function test_explicit_evaluation_state_overrides_legacy_staged_state(testCase)
   % New manifests may suppress an unavailable evaluation while retaining a
   % valid native forcing leg; only absent eval_staged fields use legacy state.
   colocation = struct();
   colocation.promice = struct('staged', true, 'eval_staged', false, ...
      'met_files', "promice/met_kanu_promice.mat", ...
      'data_files', "promice/kanu_promice.mat");

   [returned_forcing, returned_eval] = ...
      icemodel.verification.setup.colocationSourceLists(colocation);

   testCase.verifyEqual(returned_forcing, "promice");
   testCase.verifyEmpty(returned_eval);
end

function test_protocol_and_research_legs_follow_met_file_rule(testCase)
   % RetMIP/IMAU/research_site records can be eval provenance without met
   % artifacts; only staged native met files create forcing_sources.
   colocation = struct();
   colocation.retmip = struct('staged', true, ...
      'kind', 'protocol_userdata', 'data_files', "retmip/obs.mat");
   colocation.imau = struct('staged', true, ...
      'kind', 'hourly_aws', 'met_files', "met_s21_imau.mat", ...
      'forcing_ready', true);
   colocation.research_site = struct('staged', true, ...
      'kind', 'research_site_anchor');
   colocation.source_association = struct('family', "retmip", ...
      'source_id', "fa");
   colocation.nearest_noncolocated_promice = struct('nearest_anchor', "CP1");

   [returned_forcing, returned_eval] = ...
      icemodel.verification.setup.colocationSourceLists(colocation);

   testCase.verifyEqual(returned_forcing, "imau");
   expected_eval = ["retmip_protocol"; "imau_obs"; "research_site_obs"];
   testCase.verifyEqual(returned_eval, expected_eval);
end

function test_native_met_leg_without_ready_stamp_is_still_forcing_source(testCase)
   % Readiness stamps describe run-without-repair status; met files still
   % advertise selectable forcing sources when a stamp is absent.
   colocation = struct();
   colocation.promice = struct('staged', true, ...
      'met_files', "met_kanu_promice.mat", ...
      'data_files', "kanu_promice.mat");

   [returned_forcing, returned_eval] = ...
      icemodel.verification.setup.colocationSourceLists(colocation);

   testCase.verifyEqual(returned_forcing, "promice");
   testCase.verifyEqual(returned_eval, "promice_obs");
end

function test_skipped_model_legs_are_not_sources(testCase)
   % A recorded skipped leg is useful provenance, but it must not appear in
   % forcing_sources or eval_sources because no data exist to compare against.
   colocation = struct();
   colocation.sumup = struct('staged', true);
   colocation.mar = struct('staged', false, 'reason', "missing source");
   colocation.racmo = struct('staged', false, 'reason', "missing source");

   [returned_forcing, returned_eval] = ...
      icemodel.verification.setup.colocationSourceLists(colocation);

   % Skipped or data-only legs should not create forcing sources.
   expected_forcing = strings(0, 1);
   testCase.verifyEqual(returned_forcing, expected_forcing);

   % The staged SUMup observation leg is the only comparison source here.
   expected_eval = "sumup_obs";
   testCase.verifyEqual(returned_eval, expected_eval);
end

function test_rcm_leg_without_source_id_uses_product_id(testCase)
   % Source lists expose current product ids, not historical file tokens.
   colocation = struct();
   colocation.mar = struct('staged', true, ...
      'met_files', "mar/met_kanu_mar_2012_1hr.mat", ...
      'data_files', "mar/kanu_mar_2012.mat");

   [returned_forcing, returned_eval] = ...
      icemodel.verification.setup.colocationSourceLists(colocation);

   testCase.verifyEqual(returned_forcing, "mar3.11");
   testCase.verifyEqual(returned_eval, "mar3.11");
end

function test_file_bearing_legs_without_staged_flag_are_not_sources(testCase)
   % Current manifests must declare staged=true explicitly before file paths
   % become source-list entries.
   colocation = struct();
   colocation.promice = struct('met_files', "met_kanu_promice.mat");
   colocation.mar = struct('source_id', "mar3.11", ...
      'met_files', "met_kanu_mar3.11.mat", ...
      'data_files', "kanu_mar3.11.mat");

   [returned_forcing, returned_eval] = ...
      icemodel.verification.setup.colocationSourceLists(colocation);

   testCase.verifyEmpty(returned_forcing);
   testCase.verifyEmpty(returned_eval);
end

function test_unready_met_leg_remains_a_forcing_source(testCase)
   % Placeholder met artifacts are selectable forcing sources; forcing_ready
   % tells runtime code that filling or repair is required before a clean run.
   colocation = struct();
   colocation.promice = struct('staged', true, ...
      'met_files', "met_kanu_promice.mat", ...
      'data_files', "kanu_promice.mat", ...
      'forcing_ready', false, ...
      'forcing_ready_reason', "required met placeholder channel(s): ppt");

   [returned_forcing, returned_eval] = ...
      icemodel.verification.setup.colocationSourceLists(colocation);

   testCase.verifyEqual(returned_forcing, "promice");
   testCase.verifyEqual(returned_eval, "promice_obs");
end

function test_met_forcing_ready_requires_finite_required_channels(testCase)
   % A met timetable can satisfy the structural validator while still carrying
   % placeholder channels that must not be selected as runnable forcing.
   required = icemodel.forcing.helpers.metvariables();
   time = datetime(2020, 1, 1, 0, 0, 0, 'TimeZone', 'UTC') + hours((0:1)');
   met = array2timetable(ones(numel(time), numel(required)), ...
      'RowTimes', time, 'VariableNames', cellstr(required));

   met.ppt(:) = NaN;
   [ready, reason] = icemodel.verification.setup.metForcingReady(met);
   testCase.verifyFalse(ready);
   testCase.verifyTrue(contains(string(reason), "ppt"));

   met.ppt(:) = 1;
   met.ppt(1) = NaN;
   [ready, reason] = icemodel.verification.setup.metForcingReady(met);
   testCase.verifyFalse(ready);
   testCase.verifyTrue(contains(string(reason), "ppt"));

   met.ppt(:) = 1;
   [ready, reason] = icemodel.verification.setup.metForcingReady(met);
   testCase.verifyTrue(ready);
   testCase.verifyEqual(string(reason), "");

   met_missing = removevars(met, "lwd");
   [ready, reason] = icemodel.verification.setup.metForcingReady(met_missing);
   testCase.verifyFalse(ready);
   testCase.verifyTrue(contains(string(reason), "lwd"));
end

function test_met_forcing_ready_reports_contiguous_complete_windows(testCase)
   % Runnable-window discovery uses every canonical channel and splits finite
   % runs at both explicit missing values and omitted native postings.
   required = icemodel.forcing.helpers.metvariables();
   origin = datetime(2020, 1, 1, 0, 0, 0, 'TimeZone', 'UTC');
   time = origin + hours((0:5)');
   met = array2timetable(ones(numel(time), numel(required)), ...
      'RowTimes', time, 'VariableNames', cellstr(required));
   met.ppt(3) = NaN;
   met.lwd(6) = NaN;

   [ready, reason, windows] = ...
      icemodel.verification.setup.metForcingReady(met);

   testCase.verifyFalse(ready);
   testCase.verifyTrue(contains(string(reason), "ppt"));
   testCase.verifyTrue(contains(string(reason), "lwd"));
   testCase.verifyEqual(windows.start_time, time([1; 4]));
   testCase.verifyEqual(windows.end_time, time([2; 5]));
   testCase.verifyEqual(windows.sample_count, [2; 2]);

   % A missing coordinate posting also splits otherwise finite forcing.
   met = met([1, 2, 4, 5], :);
   met.ppt(:) = 1;
   met.lwd(:) = 1;
   [ready, reason, windows] = ...
      icemodel.verification.setup.metForcingReady(met);
   testCase.verifyFalse(ready);
   testCase.verifyTrue(contains(string(reason), "time coordinate"));
   testCase.verifyEqual(windows.start_time, time([1; 4]));
   testCase.verifyEqual(windows.end_time, time([2; 5]));
end

function test_met_forcing_ready_handles_empty_single_and_duplicate_times(testCase)
   % Boundary timetables stay deterministic: empty is unready, one complete
   % posting is one window, and duplicate postings are separate unready windows.
   required = icemodel.forcing.helpers.metvariables();
   origin = datetime(2020, 1, 1, 0, 0, 0, 'TimeZone', 'UTC');
   met = array2timetable(ones(1, numel(required)), ...
      'RowTimes', origin, 'VariableNames', cellstr(required));

   [ready, reason, windows] = ...
      icemodel.verification.setup.metForcingReady(met([], :));
   testCase.verifyFalse(ready);
   testCase.verifyTrue(contains(string(reason), "empty"));
   testCase.verifyEmpty(windows);

   [ready, reason, windows] = ...
      icemodel.verification.setup.metForcingReady(met);
   testCase.verifyTrue(ready);
   testCase.verifyEqual(string(reason), "");
   testCase.verifyEqual(windows.sample_count, 1);

   duplicate = array2timetable(ones(2, numel(required)), ...
      'RowTimes', [origin; origin], 'VariableNames', cellstr(required));
   [ready, reason, windows] = ...
      icemodel.verification.setup.metForcingReady(duplicate);
   testCase.verifyFalse(ready);
   testCase.verifyTrue(contains(string(reason), "time coordinate"));
   testCase.verifyEqual(windows.sample_count, [1; 1]);
end

function test_met_artifact_readiness_uses_new_and_exact_reused_bytes(testCase)
   % Readiness follows the saved path even when an exact no-overwrite call reuses
   % an earlier incomplete artifact instead of the complete request timetable.
   folder = tempname;
   mkdir(folder)
   cleanup = onCleanup(@() rmdir(folder, 's'));
   origin = datetime(2020, 1, 1, 0, 0, 0, TimeZone="UTC");
   time = origin + hours((0:23)');

   complete = readyMet(time);
   new_file = icemodel.forcing.helpers.writemet(complete, "site_new", "unit", ...
      outdir=folder, naming="window", overwrite=true);
   [ready, reason, windows] = ...
      icemodel.verification.setup.metArtifactReadiness(new_file);
   testCase.verifyTrue(ready);
   testCase.verifyEqual(string(reason), "");
   testCase.verifyEqual(numel(windows), 1);
   testCase.verifyEqual(windows.sample_count, 96);
   testCase.verifyEqual(windows.start_time, "2020-01-01T00:00:00Z");
   testCase.verifyEqual(windows.end_time, "2020-01-01T23:45:00Z");

   % A saved non-UTC display zone is converted to the same UTC instants before
   % the portable diagnostic receives its literal Z suffix.
   zoned_time = time;
   zoned_time.TimeZone = 'America/Los_Angeles';
   zoned_file = icemodel.forcing.helpers.writemet( ...
      readyMet(zoned_time), "site_zone", "unit", ...
      outdir=folder, naming="window", overwrite=true);
   [~, ~, zoned_windows] = ...
      icemodel.verification.setup.metArtifactReadiness(zoned_file);
   testCase.verifyEqual(zoned_windows.start_time, ...
      "2020-01-01T00:00:00Z");
   testCase.verifyEqual(zoned_windows.end_time, ...
      "2020-01-01T23:45:00Z");

   % Seed a separate exact artifact with an incomplete saved copy, then prove a
   % complete no-overwrite request still diagnoses that existing saved copy.
   incomplete = complete;
   incomplete.ppt(10) = NaN;
   exact_file = icemodel.forcing.helpers.writemet(incomplete, "site_exact", "unit", ...
      outdir=folder, naming="window", overwrite=true);
   returned = icemodel.forcing.helpers.writemet(complete, "site_exact", "unit", ...
      outdir=folder, naming="window", overwrite=false);
   testCase.verifyEqual(returned, exact_file);
   [ready, reason, windows] = ...
      icemodel.verification.setup.metArtifactReadiness(returned);
   testCase.verifyFalse(ready);
   testCase.verifyTrue(contains(string(reason), "ppt"));
   testCase.verifyEqual([windows.sample_count], [36 56]);

   clear cleanup
end

function test_met_artifact_readiness_uses_enclosing_reused_window(testCase)
   % A broader enclosing cache can contain a gap outside a finite narrow request;
   % manifest diagnostics must inventory the broader referenced bytes.
   folder = tempname;
   mkdir(folder)
   cleanup = onCleanup(@() rmdir(folder, 's'));
   origin = datetime(2020, 1, 1, 0, 0, 0, TimeZone="UTC");
   broad_time = origin + hours((0:47)');
   broad = readyMet(broad_time);
   broad.ppt(31) = NaN;
   broad_file = icemodel.forcing.helpers.writemet(broad, "site", "unit", ...
      outdir=folder, naming="window", overwrite=true);

   narrow = readyMet(broad_time(1:24));
   returned = icemodel.forcing.helpers.writemet(narrow, "site", "unit", ...
      outdir=folder, naming="window", overwrite=false);
   testCase.verifyEqual(returned, broad_file);

   before_bytes = binaryBytes(returned);
   before_info = dir(returned);
   [ready, reason, windows] = ...
      icemodel.verification.setup.metArtifactReadiness(returned);
   testCase.verifyFalse(ready);
   testCase.verifyTrue(contains(string(reason), "ppt"));
   testCase.verifyEqual([windows.sample_count], [120 68]);
   testCase.verifyEqual(windows(2).end_time, "2020-01-02T23:45:00Z");
   testCase.verifyEqual(binaryBytes(returned), before_bytes);
   after_info = dir(returned);
   testCase.verifyEqual(after_info.datenum, before_info.datenum);

   clear cleanup
end

function test_met_artifact_readiness_rejects_missing_and_malformed_payloads(testCase)
   % Missing paths, wrong MAT variables, structurally invalid met, and a path
   % vector all fail explicitly instead of publishing misleading readiness.
   folder = tempname;
   mkdir(folder)
   cleanup = onCleanup(@() rmdir(folder, 's'));
   missing = string(fullfile(folder, 'missing.mat'));
   testCase.verifyError( ...
      @() icemodel.verification.setup.metArtifactReadiness(missing), ...
      'icemodel:verification:metArtifactReadiness:fileMissing');

   wrong_name = fullfile(folder, 'wrong_name.mat');
   other = 1;
   save(wrong_name, 'other')
   testCase.verifyError( ...
      @() icemodel.verification.setup.metArtifactReadiness(string(wrong_name)), ...
      'icemodel:verification:metArtifactReadiness:badPayload');

   corrupt = fullfile(folder, 'corrupt.mat');
   fid = fopen(corrupt, 'w');
   cleaner = onCleanup(@() fclose(fid));
   fwrite(fid, 'not a MAT file', 'char');
   clear cleaner
   testCase.verifyError( ...
      @() icemodel.verification.setup.metArtifactReadiness(string(corrupt)), ...
      'icemodel:verification:metArtifactReadiness:loadFailed');

   bad_met = fullfile(folder, 'bad_met.mat');
   met = timetable((1:2)', 'RowTimes', ...
      datetime(2020, 1, 1, TimeZone="UTC") + hours((0:1)'), ...
      'VariableNames', {'tair'});
   save(bad_met, 'met')
   testCase.verifyError( ...
      @() icemodel.verification.setup.metArtifactReadiness(string(bad_met)), ...
      'icemodel:verification:metArtifactReadiness:badPayload');
   testCase.verifyError( ...
      @() icemodel.verification.setup.metArtifactReadiness([missing missing]), ...
      'icemodel:verification:metArtifactReadiness:scalarWindowRequired');

   clear cleanup
end

function test_write_manifest_stabilizes_forcing_complete_window_arrays(testCase)
   % Empty, singleton, and multiple nested windows remain JSON arrays after a
   % raw decode/rewrite; the unchanged singleton rewrite preserves its mtime.
   folder = tempname;
   mkdir(folder)
   cleanup = onCleanup(@() rmdir(folder, 's'));
   template = struct('start_time', "2020-01-01T00:00:00Z", ...
      'end_time', "2020-01-01T23:45:00Z", 'sample_count', 96);

   for n = 0:2
      windows = repmat(template, n, 1);
      if n == 2
         windows(2).start_time = "2020-01-02T00:00:00Z";
         windows(2).end_time = "2020-01-02T23:45:00Z";
      end
      manifest_file = fullfile(folder, sprintf('manifest_%d.json', n));
      manifest = manifestWithWindows(windows);
      icemodel.verification.setup.writeManifest(manifest_file, manifest);
      first = fileread(manifest_file);
      testCase.verifyNotEmpty(regexp(first, ...
         '"forcing_complete_windows"\s*:\s*\[', 'once'));

      % Simulate the raw-jsondecode path used by additive manifest merging.
      decoded = jsondecode(first);
      before = dir(manifest_file);
      if n == 1
         pause(1.05)
      end
      icemodel.verification.setup.writeManifest(manifest_file, decoded);
      after = dir(manifest_file);
      testCase.verifyEqual(fileread(manifest_file), first);
      if n == 1
         testCase.verifyEqual(after.datenum, before.datenum);
      end
   end

   % Already-cell input is accepted, while a nonempty scalar of another type is
   % rejected before an ambiguous manifest shape can be written.
   cell_file = fullfile(folder, 'manifest_cell.json');
   icemodel.verification.setup.writeManifest( ...
      cell_file, manifestWithWindows({template}));
   testCase.verifyNotEmpty(regexp(fileread(cell_file), ...
      '"forcing_complete_windows"\s*:\s*\[', 'once'));
   invalid = manifestWithWindows(template);
   invalid.cases.colocation.unit.forcing_complete_windows = 1;
   testCase.verifyError( ...
      @() icemodel.verification.setup.writeManifest( ...
      fullfile(folder, 'manifest_invalid.json'), invalid), ...
      'icemodel:verification:writeManifest:badForcingWindows');

   clear cleanup
end

function test_refresh_manifest_source_lists_does_not_rebuild_artifacts(testCase)
   % A source-list policy change should be repairable by rewriting manifest
   % metadata from colocation, without invoking an importer or observations
   % builder.
   folder = tempname;
   mkdir(folder)
   cleanup = onCleanup(@() rmdir(folder, 's'));
   manifest_file = fullfile(folder, 'manifest.json');

   c = struct();
   c.case_id = "kanu";
   c.forcing_sources = strings(0, 1);
   c.eval_sources = "promice_obs";
   c.colocation = struct();
   c.colocation.promice = struct('staged', false, 'eval_staged', true, ...
      'met_files', strings(1, 0), 'data_files', strings(1, 0));

   manifest = struct('dataset_family', "promice", ...
      'source_doi', "", 'source_url', "", 'source_version', "", ...
      'retrieval_date', "2026-07-09", 'cases', c, ...
      'skipped', struct([]));
   icemodel.verification.setup.writeManifest(manifest_file, manifest);

   returned = icemodel.verification.setup.refreshManifestSourceLists( ...
      manifest_file);

   testCase.verifyEmpty(returned.cases.forcing_sources);
   testCase.verifyEqual(string(returned.cases.eval_sources), "promice_obs");
   testCase.verifyTrue(isfile(manifest_file));

   text = fileread(manifest_file);
   testCase.verifyNotEmpty(regexp(text, '"forcing_sources"\s*:\s*\[', ...
      'once'));
   testCase.verifyNotEmpty(regexp(text, '"eval_sources"\s*:\s*\[', ...
      'once'));

   clear cleanup
end

function returned = leg(staged, met_file)
   %LEG Build one staged met-capable colocation leg.
   returned = struct('staged', staged, 'met_files', met_file, ...
      'data_files', "data_" + string(met_file), ...
      'forcing_ready', true);
end

function met = readyMet(time)
   %READYMET Build one finite canonical met timetable for artifact tests.
   required = icemodel.forcing.helpers.metvariables();
   met = array2timetable(ones(numel(time), numel(required)), ...
      'RowTimes', time, 'VariableNames', cellstr(required));
   met.Properties.VariableUnits = ...
      icemodel.forcing.helpers.variableUnits(required);
end

function manifest = manifestWithWindows(windows)
   %MANIFESTWITHWINDOWS Build one portable family manifest with nested windows.
   leg_metadata = struct('staged', true, ...
      'met_files', "unit/met_site_unit.mat", ...
      'forcing_ready', isempty(windows), ...
      'forcing_ready_reason', "");
   leg_metadata.forcing_complete_windows = windows;
   one_case = struct('case_id', "site", ...
      'forcing_sources', "unit", 'eval_sources', strings(0, 1), ...
      'colocation', struct('unit', leg_metadata));
   manifest = struct('dataset_family', "promice", ...
      'source_doi', "", 'source_url', "", 'source_version', "", ...
      'retrieval_date', "2026-07-12", 'cases', one_case, ...
      'skipped', struct([]));
end

function bytes = binaryBytes(filename)
   %BINARYBYTES Read exact file bytes for non-mutating helper assertions.
   fid = fopen(filename, 'r');
   cleaner = onCleanup(@() fclose(fid));
   bytes = fread(fid, Inf, '*uint8');
   clear cleaner
end
