function tests = test_promice_snow_model_ready_years
   %TEST_PROMICE_SNOW_MODEL_READY_YEARS Verify the annual handoff utility.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % Install project dependencies and build one compact synthetic family tree.
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;
   testCase.TestData.tmp = string(tempname);
   testCase.TestData.eval_root = fullfile(testCase.TestData.tmp, "eval");
   testCase.TestData.input_root = fullfile(testCase.TestData.tmp, "input");
   testCase.TestData.output_one = fullfile(testCase.TestData.tmp, "qa-one");
   testCase.TestData.output_two = fullfile(testCase.TestData.tmp, "qa-two");
   mkdir(testCase.TestData.tmp)
   writeFixtureTree(testCase.TestData.eval_root, ...
      testCase.TestData.input_root)

   % Generate once for contract assertions; the second generation is reserved
   % for the deterministic-byte test.
   testCase.TestData.report = ...
      icemodel.verification.setup.writePromiceSnowModelReadyYears( ...
      evaluation_data_root=testCase.TestData.eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      output_dir=testCase.TestData.output_one);
end

function teardownOnce(testCase)
   % Remove the isolated staged family and both generated report trees.
   if isfield(testCase.TestData, 'tmp') && isfolder(testCase.TestData.tmp)
      rmdir(testCase.TestData.tmp, 's')
   end
   if isfield(testCase.TestData, 'cleanup')
      clear testCase.TestData.cleanup
   end
end

function test_annual_scientific_contract_and_portable_outputs(testCase)
   % One fixture covers leap/non-leap counts, strict/practical evaluation,
   % daily presence, coordinate gaps, missing channels, and absent source legs.
   returned = testCase.TestData.report;
   coverage = returned.coverage;
   testCase.verifyEqual(height(coverage), 7);
   testCase.verifyEqual(height(returned.ready), 3);
   testCase.verifyEqual(nnz(returned.ready.strict_ready), 1);
   testCase.verifyTrue(returned.check.all_checks_passed);

   row = annualRecord(coverage, "MIXED", 2019);
   testCase.verifyEqual(row.expected_met_samples, 35040);
   testCase.verifyEqual(row.snow_expected_samples, 8760);
   testCase.verifyTrue(row.strict_ready);
   testCase.verifyTrue(row.promice_hybrid_practical_ready);

   row = annualRecord(coverage, "MIXED", 2020);
   testCase.verifyEqual(row.expected_met_samples, 35136);
   testCase.verifyEqual(row.snow_expected_samples, 8784);
   testCase.verifyEqual(row.tair_missing, 1);
   testCase.verifyTrue(row.practical_ready);
   testCase.verifyFalse(row.strict_ready);
   testCase.verifyFalse(row.promice_hybrid_forcing_ready);

   row = annualRecord(coverage, "FORCING_BAD", 2019);
   testCase.verifyTrue(row.rcm_time_grid_complete);
   testCase.verifyEqual(row.rcm_ppt_missing, 1);
   testCase.verifyEqual(row.rcm_required_missing_total, 2);
   testCase.verifyFalse(row.forcing_ready);

   row = annualRecord(coverage, "GRID_BAD", 2019);
   testCase.verifyFalse(row.rcm_time_grid_complete);
   testCase.verifyFalse(row.snow_time_grid_complete);
   testCase.verifyEqual(row.snow_finite_samples, 8759);
   testCase.verifyFalse(row.practical_ready);

   row = annualRecord(coverage, "DAILY_GAP", 2019);
   testCase.verifyGreaterThan(row.snow_finite_fraction, 0.95);
   testCase.verifyEqual(row.snow_finite_days, 364);
   testCase.verifyFalse(row.practical_ready);

   row = annualRecord(coverage, "PRACTICAL", 2020);
   testCase.verifyTrue(row.practical_ready);
   testCase.verifyFalse(row.strict_ready);

   row = annualRecord(coverage, "NO_SNOW", 2019);
   testCase.verifyFalse(row.snow_depth_present);
   testCase.verifyFalse(row.practical_ready);

   testCase.verifyFalse(any(coverage.site_id == "NO_DATA"));
   testCase.verifyTrue(ismember("NO_DATA", ...
      string(returned.check.sites_without_rcm_precip_leg)));

   % Generated provenance is relative and the grouped Markdown exercises both
   % a contiguous range and a practical-only site with no strict years.
   check_text = fileread(returned.files.check_json);
   summary_text = fileread(returned.files.summary_md);
   testCase.verifyFalse(contains(check_text, '/private/tmp'));
   testCase.verifySubstring(check_text, '"manifest": "promice/manifest.json"');
   testCase.verifySubstring(summary_text, '| MIXED | 2019-2020 | 2019 |');
   testCase.verifySubstring(summary_text, '| PRACTICAL | 2020 | none |');
end

function test_outputs_are_byte_deterministic_and_self_hashed(testCase)
   % Repeating the same explicit-root call must preserve every advertised hash.
   returned = icemodel.verification.setup.writePromiceSnowModelReadyYears( ...
      evaluation_data_root=testCase.TestData.eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      output_dir=testCase.TestData.output_two);
   expected = testCase.TestData.report.check;
   testCase.verifyEqual(returned.check.coverage_csv_sha256, ...
      expected.coverage_csv_sha256);
   testCase.verifyEqual(returned.check.ready_csv_sha256, ...
      expected.ready_csv_sha256);
   testCase.verifyEqual(returned.check.summary_md_sha256, ...
      expected.summary_md_sha256);
   testCase.verifyEqual( ...
      icemodel.verification.setup.fileSha256(returned.files.coverage_csv), ...
      returned.check.coverage_csv_sha256);
end

function test_missing_manifest_fails_before_writing(testCase)
   % An explicit but empty evaluation root must fail with an actionable id.
   missing_eval = fullfile(testCase.TestData.tmp, "missing-eval");
   returned = @() ...
      icemodel.verification.setup.writePromiceSnowModelReadyYears( ...
      evaluation_data_root=missing_eval, ...
      input_data_root=testCase.TestData.input_root, ...
      output_dir=fullfile(testCase.TestData.tmp, "unused"));
   testCase.verifyError(returned, ...
      'icemodel:verification:promiceReadyYears:manifestMissing');
end

%% Fixture helpers
function writeFixtureTree(eval_root, input_root)
   %WRITEFIXTURETREE Create seven annual rows spanning every acceptance branch.
   family_root = fullfile(eval_root, "promice");
   mkdir(family_root)
   mkdir(fullfile(input_root, "met", "promice"))
   mkdir(fullfile(input_root, "met", "mar3.11"))

   cases = repmat(caseTemplate(), 0, 1);
   cases(end + 1) = writeCase(family_root, input_root, "mixed", ...
      "MIXED", 2019:2020, promice_missing_year=2020, ...
      snow_missing_year=2020, snow_missing_count=1);
   cases(end + 1) = writeCase(family_root, input_root, "forcing_bad", ...
      "FORCING_BAD", 2019, mar_lwd_missing=true, ...
      mar_ppt_missing=true);
   cases(end + 1) = writeCase(family_root, input_root, "grid_bad", ...
      "GRID_BAD", 2019, drop_last_met_time=true, ...
      drop_last_snow_time=true);
   cases(end + 1) = writeCase(family_root, input_root, "daily_gap", ...
      "DAILY_GAP", 2019, snow_missing_count=24);
   cases(end + 1) = writeCase(family_root, input_root, "practical", ...
      "PRACTICAL", 2020, snow_missing_count=1);
   cases(end + 1) = writeCase(family_root, input_root, "no_snow", ...
      "NO_SNOW", 2019, write_snow=false);
   cases(end + 1) = writeCase(family_root, input_root, "no_data", ...
      "NO_DATA", 2019, write_met=false, write_snow=false);

   % The family reader expects the portable common manifest metadata.
   manifest = struct('dataset_family', 'promice', 'source_doi', '', ...
      'source_url', '', 'source_version', 'synthetic', ...
      'retrieval_date', '2026-07-12', 'cases', cases, ...
      'skipped', struct([]));
   writeJson(fullfile(family_root, "manifest.json"), manifest)
end

function c = writeCase(family_root, input_root, case_id, site_id, years, kwargs)
   %WRITECASE Save one synthetic observations bundle and optional met legs.
   arguments
      family_root (1, 1) string
      input_root (1, 1) string
      case_id (1, 1) string
      site_id (1, 1) string
      years (1, :) double
      kwargs.promice_missing_year (1, 1) double = NaN
      kwargs.snow_missing_year (1, 1) double = NaN
      kwargs.snow_missing_count (1, 1) double = 0
      kwargs.mar_lwd_missing (1, 1) logical = false
      kwargs.mar_ppt_missing (1, 1) logical = false
      kwargs.drop_last_met_time (1, 1) logical = false
      kwargs.drop_last_snow_time (1, 1) logical = false
      kwargs.write_met (1, 1) logical = true
      kwargs.write_snow (1, 1) logical = true
   end

   met_time = yearGrid(years, minutes(15));
   snow_time = yearGrid(years, hours(1));
   if kwargs.drop_last_met_time
      met_time(end) = [];
   end
   if kwargs.drop_last_snow_time
      snow_time(end) = [];
   end

   colocation = struct();
   if kwargs.write_met
      promice_met = makeMet(met_time, false);
      mar_met = makeMet(met_time, true);
      if isfinite(kwargs.promice_missing_year)
         first = find(year(promice_met.Time) == ...
            kwargs.promice_missing_year, 1);
         promice_met.tair(first) = NaN;
      end
      if kwargs.mar_lwd_missing
         mar_met.lwd(1) = NaN;
      end
      if kwargs.mar_ppt_missing
         mar_met.ppt(2) = NaN;
      end
      promice_rel = fullfile("promice", "met_" + case_id + ".mat");
      mar_rel = fullfile("mar3.11", "met_" + case_id + ".mat");
      saveMet(fullfile(input_root, "met", promice_rel), promice_met)
      saveMet(fullfile(input_root, "met", mar_rel), mar_met)
      source_window = struct( ...
         'start', sprintf('%d-01-01 00:00:00', years(1)), ...
         'end', sprintf('%d-12-31 23:00:00', years(end)));
      colocation.promice = struct('met_files', char(promice_rel), ...
         'window', source_window);
      colocation.mar = struct('met_files', char(mar_rel), ...
         'window', source_window);
   end

   data = timetable('RowTimes', snow_time);
   if kwargs.write_snow
      data.snow_depth = ones(height(data), 1);
      if kwargs.snow_missing_count > 0
         if isfinite(kwargs.snow_missing_year)
            first = find(year(data.Time) == kwargs.snow_missing_year, 1);
         else
            first = 1;
         end
         data.snow_depth(first:first + kwargs.snow_missing_count - 1) = NaN;
      end
   end
   targets = struct('format', 'timeseries', 'data', data, ...
      'metadata', struct('source', 'synthetic'));
   case_root = fullfile(family_root, case_id);
   mkdir(case_root)
   save(fullfile(case_root, "observations.mat"), 'targets')

   c = caseTemplate();
   c.case_id = char(case_id);
   c.site_id = char(site_id);
   c.surface_zone = 'ablation';
   c.period = struct('start', sprintf('%d-01-01 00:00:00', years(1)), ...
      'end', sprintf('%d-12-31 23:00:00', years(end)));
   c.evaluation_file = char(fullfile(case_id, "observations.mat"));
   c.colocation = colocation;
end

function c = caseTemplate()
   %CASETEMPLATE Keep JSON case fields homogeneous across the fixture array.
   c = struct('case_id', '', 'site_id', '', 'surface_zone', '', ...
      'period', struct('start', '', 'end', ''), 'evaluation_file', '', ...
      'colocation', struct());
end

function met = makeMet(time, include_ppt)
   %MAKEMET Create finite scalar forcing channels on the supplied coordinate.
   met = timetable('RowTimes', time);
   channels = ["tair", "swd", "lwd", "albedo", "wspd", "rh", "psfc"];
   for channel = reshape(channels, 1, [])
      met.(char(channel)) = ones(height(met), 1);
   end
   if include_ppt
      met.ppt = ones(height(met), 1);
   end
end

function time = yearGrid(years, step)
   %YEARGRID Concatenate exact UTC calendar grids for adjacent fixture years.
   time = NaT(0, 1, 'TimeZone', 'UTC');
   for y = reshape(years, 1, [])
      first = datetime(y, 1, 1, 'TimeZone', 'UTC');
      last = datetime(y + 1, 1, 1, 'TimeZone', 'UTC') - step;
      time = [time; (first:step:last)']; %#ok<AGROW>
   end
end

function saveMet(pathname, value)
   %SAVEMET Save the canonical variable name expected by the utility.
   met = value;
   save(pathname, 'met')
end

function writeJson(pathname, value)
   %WRITEJSON Save portable pretty-printed fixture metadata.
   fid = fopen(pathname, 'w');
   cleanup = onCleanup(@() fclose(fid));
   fprintf(fid, '%s\n', jsonencode(value, PrettyPrint=true));
   clear cleanup
end

function row = annualRecord(coverage, site_id, year_value)
   %ANNUALRECORD Select exactly one site/year row for concise assertions.
   keep = coverage.site_id == site_id & coverage.year == year_value;
   row = coverage(keep, :);
   if height(row) ~= 1
      error('expected exactly one annual fixture record')
   end
end
