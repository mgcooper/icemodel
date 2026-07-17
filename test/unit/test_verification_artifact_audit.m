function tests = test_verification_artifact_audit
   %TEST_VERIFICATION_ARTIFACT_AUDIT Verify reusable staged-artifact QA/QC.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   % Install project dependencies and allocate an isolated staged-data tree.
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;
   testCase.TestData.tmp = string(tempname);
   mkdir(testCase.TestData.tmp)
end

function teardown(testCase)
   % Remove generated fixture artifacts and reports after each test.
   if isfolder(testCase.TestData.tmp)
      rmdir(testCase.TestData.tmp, 's')
   end
   clear testCase.TestData.cleanup
end

function test_valid_tree_passes_writes_reports_and_remains_read_only(testCase)
   % A valid mixed-source PROMICE tree should pass while preserving explicit
   % precipitation placeholders and every staged input byte/mtime.
   [eval_root, input_root] = writeAuditTree(testCase.TestData.tmp);
   before = fileSnapshot([eval_root; input_root]);
   report_dir = fullfile(testCase.TestData.tmp, "reports");

   report = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="promice", report_dir=report_dir);

   % Print grouped diagnostics only when this nominally valid fixture regresses.
   if ~report.summary.passed
      codes = string({report.findings.code});
      for code = reshape(unique(codes, 'stable'), 1, [])
         fprintf('VALID-FIXTURE %s count=%d\n', code, nnz(codes == code));
      end
   end

   after = fileSnapshot([eval_root; input_root]);
   testCase.verifyEqual(after, before);
   testCase.verifyTrue(report.summary.passed);
   testCase.verifyEqual(report.summary.error_count, 0);
   testCase.verifyEqual(report.summary.blocker_count, 0);
   testCase.verifyTrue(any(string({report.findings.severity}) == "placeholder"));
   testCase.verifyTrue(isfile(report.report_files.json));
   testCase.verifyTrue(isfile(report.report_files.markdown));
   decoded = jsondecode(fileread(report.report_files.json));
   testCase.verifyTrue(decoded.summary.passed);
   testCase.verifyTrue(isfield(decoded.artifacts, 'artifact_sha256'));
   testCase.verifyTrue(isfield(decoded.artifacts, 'artifact_size_bytes'));
   testCase.verifySubstring(fileread(report.report_files.markdown), ...
      "# Verification artifact QA/QC");

   % MAR magnitude diagnostics remain informational while UTC-day reconstruction
   % is checked directly.
   runoff = channelRecord(report.channels, "mar3.11", "runoff");
   testCase.verifyEqual(runoff.nonconstant_utc_day_count, 0);
   testCase.verifyGreaterThan(runoff.complete_utc_day_count, 0);
   testCase.verifyEqual(runoff.partial_utc_day_count, 2);
   testCase.verifyTrue(isfinite(runoff.p99));
   subl = channelRecord(report.channels, "mar3.11", "subl");
   testCase.verifyLessThan(subl.minimum, 0);
   testCase.verifyGreaterThan(subl.maximum, 0);
end

function test_audits_mar_profile_model_output_contract(testCase)
   % Optional MAR RO1 sidecars pass with complete schema/provenance and fail
   % when the public density-unit contract is altered.
   [eval_root, input_root, paths] = writeAuditTree(testCase.TestData.tmp);
   profile_file = attachMarProfileFixture(input_root, paths.manifest);

   report = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="promice", report_dir="");
   profile_record = report.artifacts( ...
      string({report.artifacts.kind}) == "model_output");
   codes = string({report.findings.code});
   testCase.verifyEqual(string(profile_record.status), "passed")
   testCase.verifyFalse(any(startsWith(codes, "profile_")))

   loaded = load(profile_file, 'reference');
   density_index = strcmp( ...
      loaded.reference.data.density.Properties.VariableNames, 'density');
   loaded.reference.data.density.Properties.VariableUnits{density_index} = 'kg';
   reference = loaded.reference;
   save(profile_file, 'reference')
   invalid = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="promice", report_dir="");
   testCase.verifyTrue(any(string({invalid.findings.code}) == "profile_units"))
end

function test_retMIP_netcdf_and_status_only_manifest_contract(testCase)
   % Native RetMIP protocol NetCDF files are inspected as NetCDF, while the
   % native_met readiness leg remains metadata rather than a missing artifact.
   [eval_root, input_root, protocol_file] = ...
      writeRetmipAuditTree(testCase.TestData.tmp);

   report = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="retmip", report_dir="");
   codes = string({report.findings.code});
   protocol = report.artifacts( ...
      string({report.artifacts.path}) == protocol_file);
   testCase.verifyEqual(string(protocol.payload), "netcdf");
   testCase.verifyEqual(protocol.n_channels, 2);
   testCase.verifyEqual(protocol.n_samples, 3);
   testCase.verifyEqual(string(protocol.status), "passed");
   testCase.verifyFalse(any(ismember(codes, ["artifact_read_error", ...
      "missing_artifact", "staged_leg_without_artifact", ...
      "eval_source_without_artifact"])));

   % Native participant time conventions remain source-faithful: t aliases,
   % absent coordinate variables/units, and noncanonical values are not
   % relabeled by the audit.
   writeRetmipProtocolNetcdf(protocol_file, "missing_time");
   missing_time = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="retmip", report_dir="");
   testCase.verifyFalse(any(string({missing_time.findings.code}) ...
      == "artifact_read_error"));

   writeRetmipProtocolNetcdf(protocol_file, "missing_units");
   missing_units = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="retmip", report_dir="");
   testCase.verifyFalse(any(string({missing_units.findings.code}) ...
      == "artifact_read_error"));

   writeRetmipProtocolNetcdf(protocol_file, "nonmonotonic");
   nonmonotonic = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="retmip", report_dir="");
   testCase.verifyFalse(any(string({nonmonotonic.findings.code}) ...
      == "nonmonotonic_time"));

   writeRetmipProtocolNetcdf(protocol_file, "nonfinite");
   nonfinite = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="retmip", report_dir="");
   testCase.verifyFalse(any(string({nonfinite.findings.code}) ...
      == "missing_time"));

   % A NetCDF container still has to carry the expected protocol payload.
   writeRetmipProtocolNetcdf(protocol_file, "missing_variables");
   missing_variables = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="retmip", report_dir="");
   testCase.verifyTrue(any(string({missing_variables.findings.code}) ...
      == "netcdf_protocol_variables"));
end

function test_audits_transient_albedo_provenance_and_derived_mask(testCase)
   % The audit recomputes transient episodes from preserved raw radiation and
   % rejects either finite derived values or stale compact provenance.
   [eval_root, input_root, ~] = ...
      writeRetmipAuditTree(testCase.TestData.tmp);
   [data_file, flagged] = attachAlbedoTransientFixture( ...
      eval_root, input_root);

   valid = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="retmip", report_dir="");
   valid_codes = string({valid.findings.code});
   testCase.verifyFalse(any(startsWith( ...
      valid_codes, "albedo_transient_qc_")));

   % One leaked derived channel is sufficient to make the artifact invalid;
   % the raw swd/swu values remain unchanged.
   loaded = load(data_file, 'Data');
   loaded.Data.albedo(flagged) = 0.4;
   saveData(data_file, loaded.Data)
   leaked = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="retmip", report_dir="");
   testCase.verifyTrue(any(string({leaked.findings.code}) ...
      == "albedo_transient_qc_mask"));

   % Correct masking cannot compensate for a stale count stamped into both
   % timetable and top-level artifact metadata.
   loaded.Data.albedo(flagged) = NaN;
   loaded.Data.Properties.UserData.albedo_transient_qc.episode_day_count = ...
      loaded.Data.Properties.UserData.albedo_transient_qc.episode_day_count + 1;
   saveData(data_file, loaded.Data)
   stale = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="retmip", report_dir="");
   testCase.verifyTrue(any(string({stale.findings.code}) ...
      == "albedo_transient_qc_contract"));
end

function test_profile_event_times_units_and_ranges(testCase)
   % SUMup-style profile/event rows may repeat or arrive out of order, but
   % uncertainty units, usable event dates, and source-unit ranges stay strict.
   [eval_root, input_root, observations_file] = ...
      writeProfileAuditTree(testCase.TestData.tmp);

   report = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="sumup", report_dir="");
   codes = string({report.findings.code});
   testCase.verifyFalse(any(ismember(codes, ["nonmonotonic_time", ...
      "missing_time", "missing_unit", "wrong_unit", "physical_range"])));

   % Generic uncertainty metadata is interpreted through its parent channel;
   % blank and mismatched units remain real artifact defects.
   loaded = load(observations_file, 'targets');
   loaded.targets.data.density.Properties.VariableUnits{1} = '';
   targets = loaded.targets;
   save(observations_file, 'targets')
   blank = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="sumup", report_dir="");
   keep = string({blank.findings.code}) == "missing_unit" ...
      & string({blank.findings.channel}) == "error";
   testCase.verifyTrue(any(keep));
   testCase.verifySubstring(string(blank.findings(find(keep, 1)).message), ...
      "kg m-3");

   loaded.targets.data.density.Properties.VariableUnits{1} = 'm';
   targets = loaded.targets;
   save(observations_file, 'targets')
   wrong_error = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="sumup", report_dir="");
   keep = string({wrong_error.findings.code}) == "wrong_unit" ...
      & string({wrong_error.findings.channel}) == "error";
   testCase.verifyTrue(any(keep));

   % degC is a valid source unit for observational profile temperatures, but an
   % unrelated unit or physically impossible degC value is not.
   loaded.targets.data.density.Properties.VariableUnits{1} = 'kg m-3';
   loaded.targets.data.subsurface_temperature.Properties.VariableUnits{2} = 'm';
   targets = loaded.targets;
   save(observations_file, 'targets')
   wrong_temperature = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="sumup", report_dir="");
   keep = string({wrong_temperature.findings.code}) == "wrong_unit" ...
      & string({wrong_temperature.findings.channel}) ...
      == "subsurface_temperature";
   testCase.verifyTrue(any(keep));

   loaded.targets.data.subsurface_temperature.Properties.VariableUnits{2} = ...
      'degC';
   loaded.targets.data.subsurface_temperature.subsurface_temperature(1) = 100;
   targets = loaded.targets;
   save(observations_file, 'targets')
   invalid_temperature = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="sumup", report_dir="");
   keep = string({invalid_temperature.findings.code}) == "physical_range" ...
      & string({invalid_temperature.findings.channel}) ...
      == "subsurface_temperature";
   testCase.verifyTrue(any(keep));

   % An interval row with neither a date nor fallback year remains invalid.
   loaded.targets.data.subsurface_temperature.subsurface_temperature(1) = -20;
   loaded.targets.data.smb.start_year(1) = NaN;
   targets = loaded.targets;
   save(observations_file, 'targets')
   invalid_event = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="sumup", report_dir="");
   testCase.verifyTrue(any(string({invalid_event.findings.code}) ...
      == "missing_time"));
end

function test_model_profiles_use_observation_not_forcing_window(testCase)
   % Exact-date MAR profiles follow observation dates even when the reusable
   % forcing leg begins later than an otherwise valid profile survey.
   [eval_root, input_root, paths] = writeAuditTree(testCase.TestData.tmp);
   attachMarProfileFixture(input_root, paths.manifest);
   manifest = jsondecode(fileread(paths.manifest));
   manifest.cases.colocation.mar.window.start = '2000-06-01 00:00:00';
   writeJson(paths.manifest, manifest)

   report = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="promice", report_dir="");
   profile = string({report.artifacts.kind}) == "model_output";
   testCase.verifyTrue(any(profile));
   testCase.verifyFalse(any(string({report.findings.code}) ...
      == "period_outside_manifest"));
end

function test_records_missing_runs_and_immutable_artifact_identity(testCase)
   % Missing-run metrics retain the existing missing total while distinguishing
   % separate outages. File identity must describe the exact audited bytes.
   [eval_root, input_root, paths] = writeAuditTree(testCase.TestData.tmp);
   loaded = load(paths.promice_data, 'Data');
   loaded.Data.tair = [NaN; 260; NaN; NaN];
   saveData(paths.promice_data, loaded.Data)

   report = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="promice", report_dir="");

   % The three missing samples form two independent hourly runs, while a
   % complete met channel explicitly reports no runs.
   channels = report.channels;
   data_keep = string({channels.source}) == "promice" ...
      & string({channels.kind}) == "userdata" ...
      & string({channels.channel}) == "tair";
   data_tair = channels(find(data_keep, 1));
   testCase.verifyEqual(data_tair.missing_count, 3);
   testCase.verifyEqual(data_tair.missing_run_count, 2);
   testCase.verifyEqual(data_tair.longest_missing_run_samples, 2);
   met_keep = string({channels.source}) == "promice" ...
      & string({channels.kind}) == "met" ...
      & string({channels.channel}) == "tair";
   met_tair = channels(find(met_keep, 1));
   testCase.verifyEqual(met_tair.missing_run_count, 0);
   testCase.verifyEqual(met_tair.longest_missing_run_samples, 0);

   % Every present ledger entry carries the shared SHA-256 and exact byte size.
   artifacts = report.artifacts;
   present = artifacts([artifacts.exists]);
   testCase.verifyTrue(all(strlength(string({present.artifact_sha256})) == 64));
   testCase.verifyTrue(all([present.artifact_size_bytes] > 0));
   path_keep = string({artifacts.path}) == string(paths.promice_data);
   audited = artifacts(find(path_keep, 1));
   current_info = dir(paths.promice_data);
   testCase.verifyEqual(audited.artifact_size_bytes, current_info.bytes);
   testCase.verifyEqual(audited.artifact_sha256, ...
      icemodel.verification.setup.fileSha256(paths.promice_data));

   % A same-size post-audit byte mutation proves SHA-256 catches replacements
   % that byte-size checks alone cannot distinguish. The immutable report record
   % continues to describe only the audited file.
   fid = fopen(paths.promice_data, 'r+');
   testCase.assertGreaterThanOrEqual(fid, 0);
   cleanup = onCleanup(@() fclose(fid));
   first_byte = fread(fid, 1, '*uint8');
   testCase.assertEqual(numel(first_byte), 1);
   fseek(fid, 0, 'bof');
   fwrite(fid, bitxor(first_byte, uint8(1)), 'uint8');
   clear cleanup
   changed_info = dir(paths.promice_data);
   testCase.verifyEqual(changed_info.bytes, audited.artifact_size_bytes);
   testCase.verifyNotEqual( ...
      icemodel.verification.setup.fileSha256(paths.promice_data), ...
      audited.artifact_sha256);
end

function test_mar_15m_met_uses_derived_not_native_daily_contract(testCase)
   % A manifest-referenced MAR met artifact keeps the native-ledger provenance
   % copied from Data, but the audit must judge its 15-minute payload only by
   % generic schema/time, metadata-sync, and zero-order-hold contracts.
   [eval_root, input_root, paths] = writeAuditTree(testCase.TestData.tmp);
   met_file = attachMarMetFixture(input_root, paths.manifest, paths.mar_data);

   report = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="promice", report_dir="");
   codes = string({report.findings.code});
   native_codes = ["mar_qc_count_mismatch", "mar_qc_day_ledger_invalid", ...
      "mar_daily_constraint_inconsistent", "mar_melt_daily_validation", ...
      "mar_diagnostic_metadata_missing"];
   testCase.verifyTrue(report.summary.passed, strjoin(codes, newline));
   testCase.verifyFalse(any(ismember(codes, native_codes)));

   % Independent corruptions prove that narrowing the MAR-native checks does
   % not weaken the contracts shared by every 15-minute met artifact.
   loaded = load(met_file, 'met', 'artifact_metadata');
   runoff_index = find(string(loaded.met.Properties.VariableNames) == "runoff", 1);
   loaded.met.Properties.VariableUnits{runoff_index} = 'invalid';
   loaded.met.Time(2) = loaded.met.Time(1);
   metadata = loaded.met.Properties.UserData;
   metadata.met_resample_expected_missing_counts.runoff = 1;
   loaded.met.Properties.UserData = metadata;
   met = loaded.met;
   artifact_metadata = metadata;
   artifact_metadata.sample_method = 'bilinear';
   save(met_file, 'met', 'artifact_metadata')

   bad = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="promice", report_dir="");
   bad_codes = string({bad.findings.code});
   testCase.verifyTrue(all(ismember(["wrong_unit", "nonmonotonic_time", ...
      "met_gap_interpolation", "artifact_metadata_mismatch"], bad_codes)), ...
      strjoin(bad_codes, newline));

   % The derived/native distinction is not permission to ignore manifest refs.
   delete(met_file)
   missing = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="promice", report_dir="");
   testCase.verifyTrue(any(string({missing.findings.code}) == "missing_artifact"));
end

function test_detects_schema_time_range_unit_period_and_metadata_errors(testCase)
   % One malformed met file should surface independent schema, time, range,
   % period, and metadata-parity failures in a single non-throwing pass.
   [eval_root, input_root, paths] = writeAuditTree(testCase.TestData.tmp);
   loaded = load(paths.promice_met, 'met', 'artifact_metadata');
   loaded.met.Time(2) = loaded.met.Time(1);
   loaded.met.albedo(1) = 1.5;
   swd = string(loaded.met.Properties.VariableNames) == "swd";
   loaded.met.Properties.VariableUnits{swd} = 'bad-unit';
   loaded.met.Properties.VariableNames{ ...
      string(loaded.met.Properties.VariableNames) == "tair"} = 'bad_tair';
   artifact_metadata = loaded.artifact_metadata;
   artifact_metadata.stale_top_level = true;
   met = loaded.met;
   save(paths.promice_met, 'met', 'artifact_metadata')

   manifest = jsondecode(fileread(paths.manifest));
   manifest.cases.colocation.promice.window.start = '1999-12-31 23:00:00';
   writeJson(paths.manifest, manifest)

   report = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="promice");
   codes = string({report.findings.code});

   testCase.verifyFalse(report.summary.passed);
   verifyCodes(testCase, codes, ["met_contract", "noncanonical_channel", ...
      "wrong_unit", "physical_range", "nonmonotonic_time", ...
      "period_coverage_drift", "artifact_metadata_mismatch"]);
end

function test_detects_filled_met_gap_and_missing_resample_provenance(testCase)
   % Source-derived missing-count proof catches a finite derived-met bridge.
   [eval_root, input_root, paths] = writeAuditTree(testCase.TestData.tmp);
   loaded = load(paths.promice_met, 'met', 'artifact_metadata');
   loaded.met.Properties.UserData.met_resample_expected_missing_counts.tair = 4;
   loaded.artifact_metadata = loaded.met.Properties.UserData;
   met = loaded.met;
   artifact_metadata = loaded.artifact_metadata;
   save(paths.promice_met, 'met', 'artifact_metadata')

   report = icemodel.verification.auditArtifacts( ...
      families="promice", evaluation_data_root=eval_root, ...
      input_data_root=input_root);
   testCase.verifyTrue(any(string({report.findings.code}) ...
      == "met_gap_interpolation"));

   % Removing the durable source-gap proof is independently actionable.
   metadata = met.Properties.UserData;
   metadata = rmfield(metadata, {'met_resample_policy', ...
      'met_resample_expected_missing_counts'});
   met.Properties.UserData = metadata;
   artifact_metadata = metadata;
   save(paths.promice_met, 'met', 'artifact_metadata')
   report = icemodel.verification.auditArtifacts( ...
      families="promice", evaluation_data_root=eval_root, ...
      input_data_root=input_root);
   testCase.verifyTrue(any(string({report.findings.code}) ...
      == "met_resample_gap_metadata_missing"));
end

function test_detects_merra_mar_modis_promice_and_missing_artifact(testCase)
   % Source-specific defects should retain distinct codes so downstream reports
   % can route repairs without parsing prose.
   [eval_root, input_root, paths] = writeAuditTree(testCase.TestData.tmp);

   % MERRA's numeric channels are invalid when durable metadata still declares
   % the native positive-upward orientation.
   loaded = load(paths.merra_data, 'Data', 'artifact_metadata');
   loaded.Data.Properties.UserData.merra_flux_sign_convention = 'positive_upward';
   loaded.artifact_metadata = loaded.Data.Properties.UserData;
   Data = loaded.Data;
   artifact_metadata = loaded.artifact_metadata;
   save(paths.merra_data, 'Data', 'artifact_metadata')

   % Contradict MODIS year metadata while retaining finite year-2000 values.
   loaded = load(paths.merra_data, 'Data', 'artifact_metadata');
   loaded.Data.Properties.UserData.modis_coverage_years = 1999;
   loaded.artifact_metadata = loaded.Data.Properties.UserData;
   Data = loaded.Data;
   artifact_metadata = loaded.artifact_metadata;
   save(paths.merra_data, 'Data', 'artifact_metadata')

   % One altered hour violates native-daily constancy; independent metadata
   % mutations exercise source policy, channel, and day-count consistency.
   loaded = load(paths.mar_data, 'Data', 'artifact_metadata');
   loaded.Data.runoff(3) = loaded.Data.runoff(3) + 0.01;
   loaded.Data.Properties.UserData.mar_qc_hourly_distribution = 'wrong';
   loaded.Data.Properties.UserData.mar_qc_channels = "runoff";
   loaded.Data.Properties.UserData.mar_qc_complete_utc_day_count = ...
      loaded.Data.Properties.UserData.mar_qc_complete_utc_day_count + 1;
   loaded.Data.Properties.UserData.mar_snowd_masked_years = 2000;
   Data = loaded.Data;
   artifact_metadata = loaded.Data.Properties.UserData;
   save(paths.mar_data, 'Data', 'artifact_metadata')

   % PROMICE must never promote liquid gauge data into total precipitation.
   loaded = load(paths.promice_met, 'met', 'artifact_metadata');
   loaded.met.ppt(1) = 1e-8;
   met = loaded.met;
   artifact_metadata = loaded.artifact_metadata;
   save(paths.promice_met, 'met', 'artifact_metadata')

   % A deleted target proves exact manifest/file existence checking.
   delete(paths.observations)

   report = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="promice");
   codes = string({report.findings.code});

   verifyCodes(testCase, codes, ["merra_sign_not_canonical", ...
      "mar_daily_constraint_inconsistent", "mar_qc_metadata_invalid", ...
      "mar_qc_channels_mismatch", "mar_qc_count_mismatch", ...
      "mar_snowd_mask_invalid", ...
      "modis_coverage_outside_artifact", ...
      "modis_coverage_data_mismatch", "promice_invented_precip", ...
      "missing_artifact"]);
end

function test_audits_mar_optional_diagnostic_contract(testCase)
   % Source-light QA rejects wrong source semantics, broken daily support,
   % stale signed-RZ statistics, unconverted magnitudes, and an ME/MEH mismatch.
   [eval_root, input_root, paths] = writeAuditTree(testCase.TestData.tmp);
   loaded = load(paths.mar_data, 'Data', 'artifact_metadata');
   Data = loaded.Data;
   Data.Properties.UserData.mar_diagnostic_subl_sign = 'wrong';
   Data.subl_evap(3) = Data.subl_evap(3) + 0.001;
   Data.refreeze_deposition(4) = -0.1;
   Data.melt(2) = Data.melt(2) + 0.01;
   artifact_metadata = Data.Properties.UserData;
   save(paths.mar_data, 'Data', 'artifact_metadata')

   report = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="promice");
   codes = string({report.findings.code});

   verifyCodes(testCase, codes, ["mar_diagnostic_metadata_invalid", ...
      "mar_diagnostic_range", "mar_refreeze_metadata_mismatch", ...
      "mar_daily_diagnostic_support", "mar_melt_daily_validation"]);
end

function test_accepts_signed_native_mar_refreeze_with_matching_metadata(testCase)
   % Source-exact strict and material RZ negatives are valid when both durable
   % statistic sets and constant native daily support agree with userdata.
   [eval_root, input_root, paths] = writeAuditTree(testCase.TestData.tmp);
   loaded = load(paths.mar_data, 'Data', 'artifact_metadata');
   Data = loaded.Data;
   days = dateshift(Data.Time, 'start', 'day');
   selected = unique(days, 'stable');
   minima = [-1.46520137786865e-5, ...
      -7.63545682032903e-6, -5.54515173037847e-6];
   for k = 1:numel(minima)
      Data.refreeze_deposition(days == selected(k)) = minima(k);
   end
   Data.refreeze_deposition(days == selected(4)) = -1e-10;
   artifact_metadata = icemodel.forcing.helpers.marRefreezeMetadata( ...
      Data, Data.Properties.UserData);
   Data.Properties.UserData = artifact_metadata;
   save(paths.mar_data, 'Data', 'artifact_metadata')

   report = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="promice");
   codes = string({report.findings.code});

   testCase.verifyFalse(any(ismember(codes, [ ...
      "mar_refreeze_metadata_mismatch", "mar_diagnostic_range", ...
      "mar_daily_diagnostic_support"])));
end

function test_rejects_tampered_signed_mar_refreeze_metadata(testCase)
   % Policy, negative-day count, and minimum are independently durable rather
   % than advisory labels that can drift away from the saved userdata values.
   [eval_root, input_root, paths] = writeAuditTree(testCase.TestData.tmp);
   loaded = load(paths.mar_data, 'Data', 'artifact_metadata');
   Data = loaded.Data;
   days = dateshift(Data.Time, 'start', 'day');
   first_day = min(days);
   Data.refreeze_deposition(days == first_day) = -1e-4;
   metadata = icemodel.forcing.helpers.marRefreezeMetadata( ...
      Data, Data.Properties.UserData);
   metadata.mar_diagnostic_refreeze_deposition_sign = 'wrong';
   metadata.mar_diagnostic_refreeze_negative_day_count = 2;
   metadata.mar_diagnostic_refreeze_negative_minimum_mwe_h = -2e-4;
   metadata.mar_diagnostic_refreeze_material_negative_day_count = 2;
   metadata.mar_diagnostic_refreeze_material_negative_minimum_mwe_h = -3e-4;
   Data.Properties.UserData = metadata;
   artifact_metadata = metadata;
   save(paths.mar_data, 'Data', 'artifact_metadata')

   report = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="promice");
   codes = string({report.findings.code});

   verifyCodes(testCase, codes, ["mar_diagnostic_metadata_invalid", ...
      "mar_refreeze_metadata_mismatch"]);
end

function test_accepts_reduced_mar_optional_diagnostic_schema(testCase)
   % Removing every optional native product remains a valid, explicit reduced
   % schema and does not turn the MAR model-met leg into a readiness failure.
   [eval_root, input_root, paths] = writeAuditTree(testCase.TestData.tmp);
   loaded = load(paths.mar_data, 'Data', 'artifact_metadata');
   Data = removevars(loaded.Data, ...
      {'subl', 'subl_evap', 'refreeze_deposition'});
   metadata = loaded.Data.Properties.UserData;
   fields = string(fieldnames(metadata));
   diagnostic_fields = fields(startsWith(fields, "mar_diagnostic_"));
   metadata = rmfield(metadata, cellstr(diagnostic_fields));
   metadata = icemodel.forcing.helpers.marDiagnosticMetadata( ...
      Data, [], metadata, sector=1);
   Data.Properties.UserData = metadata;
   artifact_metadata = metadata;
   save(paths.mar_data, 'Data', 'artifact_metadata')

   report = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="promice");
   codes = string({report.findings.code});
   testCase.verifyFalse(any(startsWith(codes, "mar_diagnostic_") ...
      | codes == "mar_melt_daily_validation"));
   testCase.verifyTrue(report.summary.passed);
end

function test_rejects_missing_or_inconsistent_mar_diagnostic_inventory(testCase)
   % Missing provenance and a forged source/channel inventory have distinct
   % repair codes so neither can masquerade as an optional reduced source.
   [eval_root, input_root, paths] = writeAuditTree(testCase.TestData.tmp);
   loaded = load(paths.mar_data, 'Data', 'artifact_metadata');
   original = loaded.Data.Properties.UserData;
   Data = loaded.Data;
   metadata = rmfield(original, 'mar_diagnostic_method');
   Data.Properties.UserData = metadata;
   artifact_metadata = metadata;
   save(paths.mar_data, 'Data', 'artifact_metadata')
   report = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="promice");
   testCase.verifyTrue(any(string({report.findings.code}) ...
      == "mar_diagnostic_metadata_missing"));

   metadata = original;
   metadata.mar_diagnostic_channels = "subl";
   Data.Properties.UserData = metadata;
   artifact_metadata = metadata;
   save(paths.mar_data, 'Data', 'artifact_metadata')
   report = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="promice");
   testCase.verifyTrue(any(string({report.findings.code}) ...
      == "mar_diagnostic_channels_mismatch"));
end

function test_audits_promice_tice10m_source_mask_and_provenance(testCase)
   % The reusable audit reports safely masked thermistor failures and rejects
   % either a finite flagged target or incomplete durable QC provenance.
   [eval_root, input_root, paths] = writeAuditTree(testCase.TestData.tmp);
   loaded = load(paths.promice_data, 'Data', 'artifact_metadata');
   Data = loaded.Data;
   source = [260; 260.1; 267; 267.1];
   flag = [0; 1; 1; 0];
   target = source;
   target(flag > 0) = NaN;
   Data.tice10m = target;
   Data.tice10m_source = source;
   Data.tice10m_qc_flag = flag;
   Data = icemodel.forcing.helpers.stampMetadata(Data);
   metadata = Data.Properties.UserData;
   metadata.tice10m_qc_status = "applied";
   metadata.tice10m_qc_method = ...
      "mask_gt_1K_hourly_endpoints_and_large_isolated_sensor_epochs";
   metadata.tice10m_qc_source_variable = "t_i_10m";
   metadata.tice10m_qc_source_channel = "tice10m_source";
   metadata.tice10m_qc_jump_threshold_K = 1;
   metadata.tice10m_qc_persistent_jump_threshold_K = 4;
   metadata.tice10m_qc_other_sensor_median_threshold_K = 0.25;
   metadata.tice10m_qc_target_depth_tolerance_m = 2;
   metadata.tice10m_qc_recovery_window_hours = 24;
   metadata.tice10m_qc_depth_reset_threshold_m = 0.5;
   metadata.tice10m_qc_flag_codes = struct( ...
      'accepted', 0, 'failed', 1, 'unreviewed', 2, ...
      'persistent_unreviewed', 3);
   metadata.tice10m_qc_flagged_sample_count = 2;
   metadata.tice10m_qc_failed_sample_count = 2;
   metadata.tice10m_qc_unreviewed_sample_count = 0;
   metadata.tice10m_qc_persistent_sample_count = 0;
   metadata.tice10m_qc_basis = "synthetic reviewed discontinuity";
   Data.Properties.UserData = metadata;
   artifact_metadata = metadata;
   save(paths.promice_data, 'Data', 'artifact_metadata')

   report = icemodel.verification.auditArtifacts( ...
      families="promice", evaluation_data_root=eval_root, ...
      input_data_root=input_root);
   codes = string({report.findings.code});
   testCase.verifyTrue(any(codes == ...
      "promice_tice10m_failed_samples_masked"));
   testCase.verifyFalse(any(codes == "promice_tice10m_qc_mask_invalid"));

   % Recompute the temporal invariant from the preserved source rather than
   % trusting a self-consistent but falsely all-accepted flag/count record.
   Data.tice10m = source;
   Data.tice10m_qc_flag = zeros(size(flag));
   metadata.tice10m_qc_flagged_sample_count = 0;
   metadata.tice10m_qc_failed_sample_count = 0;
   Data.Properties.UserData = metadata;
   artifact_metadata = metadata;
   save(paths.promice_data, 'Data', 'artifact_metadata')
   report = icemodel.verification.auditArtifacts( ...
      families="promice", evaluation_data_root=eval_root, ...
      input_data_root=input_root);
   testCase.verifyTrue(any(string({report.findings.code}) ...
      == "promice_tice10m_discontinuity_unmasked"));

   % Cadence-independent source-range QC must reject a forged accepted one-row
   % artifact; there is no temporal pair from which the jump audit could infer it.
   full_Data = Data;
   Data = Data(1, :);
   Data.tice10m_source = icemodel.physicalConstant('Tf') + 2;
   Data.tice10m = Data.tice10m_source;
   Data.tice10m_qc_flag = 0;
   Data.Properties.UserData = metadata;
   artifact_metadata = metadata;
   save(paths.promice_data, 'Data', 'artifact_metadata')
   report = icemodel.verification.auditArtifacts( ...
      families="promice", evaluation_data_root=eval_root, ...
      input_data_root=input_root);
   testCase.verifyTrue(any(string({report.findings.code}) ...
      == "promice_tice10m_source_range_unmasked"));
   Data = full_Data;

   % Code 3 is a separately reported unresolved epoch while remaining safely
   % masked under the same canonical source/target contract.
   Data.tice10m = target;
   Data.tice10m_qc_flag = [0; 3; 3; 0];
   metadata.tice10m_qc_flagged_sample_count = 2;
   metadata.tice10m_qc_failed_sample_count = 0;
   metadata.tice10m_qc_unreviewed_sample_count = 2;
   metadata.tice10m_qc_persistent_sample_count = 2;
   Data.Properties.UserData = metadata;
   artifact_metadata = metadata;
   save(paths.promice_data, 'Data', 'artifact_metadata')
   report = icemodel.verification.auditArtifacts( ...
      families="promice", evaluation_data_root=eval_root, ...
      input_data_root=input_root);
   testCase.verifyTrue(any(string({report.findings.code}) ...
      == "promice_tice10m_persistent_epoch_masked"));

   % Count provenance is part of the durable contract and cannot silently drift
   % from the staged flag after an additive refresh or daily aggregation.
   metadata.tice10m_qc_flagged_sample_count = 99;
   Data.Properties.UserData = metadata;
   artifact_metadata = metadata;
   save(paths.promice_data, 'Data', 'artifact_metadata')
   report = icemodel.verification.auditArtifacts( ...
      families="promice", evaluation_data_root=eval_root, ...
      input_data_root=input_root);
   testCase.verifyTrue(any(string({report.findings.code}) ...
      == "promice_tice10m_qc_count_mismatch"));
   metadata.tice10m_qc_flagged_sample_count = 2;
   Data.Properties.UserData = metadata;
   artifact_metadata = metadata;

   % A finite canonical value at a flagged endpoint violates the core safety
   % contract even though the source and metadata still claim successful QC.
   Data.tice10m(2) = Data.tice10m_source(2);
   save(paths.promice_data, 'Data', 'artifact_metadata')
   report = icemodel.verification.auditArtifacts( ...
      families="promice", evaluation_data_root=eval_root, ...
      input_data_root=input_root);
   testCase.verifyTrue(any(string({report.findings.code}) ...
      == "promice_tice10m_qc_mask_invalid"));

   % Removing one required field proves stale pre-contract artifacts cannot pass
   % merely by carrying channels with plausible names and numeric values.
   metadata = rmfield(metadata, 'tice10m_qc_method');
   Data.Properties.UserData = metadata;
   artifact_metadata = metadata;
   save(paths.promice_data, 'Data', 'artifact_metadata')
   report = icemodel.verification.auditArtifacts( ...
      families="promice", evaluation_data_root=eval_root, ...
      input_data_root=input_root);
   testCase.verifyTrue(any(string({report.findings.code}) ...
      == "promice_tice10m_qc_contract_missing"));
end

function test_audits_hourly_promice_shortwave_darkness_only(testCase)
   % Missing deep-night zeros are invalid only in hourly native userdata; the
   % same timestamp in 15-minute met remains governed by resampling provenance.
   [eval_root, input_root, paths] = writeAuditTree(testCase.TestData.tmp);
   loaded = load(paths.promice_data, 'Data', 'artifact_metadata');
   Data = loaded.Data;
   Data.Properties.DimensionNames{1} = 'Time';
   Data.swd = [0; NaN; 0; 0];
   Data.swu = [0; NaN; 0; 0];
   Data = icemodel.forcing.helpers.stampMetadata(Data);
   metadata = Data.Properties.UserData;
   metadata.lat = 70;
   metadata.lon = -40;
   Data.Properties.UserData = metadata;
   artifact_metadata = metadata;
   hourly_Data = Data;
   save(paths.promice_data, 'Data', 'artifact_metadata')

   report = icemodel.verification.auditArtifacts( ...
      families="promice", evaluation_data_root=eval_root, ...
      input_data_root=input_root);
   codes = string({report.findings.code});
   gaps = report.findings(codes == "promice_shortwave_darkness_gap");
   testCase.verifyEqual(sort(string({gaps.channel})), ["swd", "swu"]);

   % A surgical all-NaN window still needs deep-night zeros when whole-file
   % metadata proves that each source channel has observations elsewhere.
   Data.swd(:) = NaN;
   Data.swu(:) = NaN;
   metadata.swd_source_file_observations_present = true;
   metadata.swu_source_file_observations_present = true;
   Data.Properties.UserData = metadata;
   artifact_metadata = metadata;
   save(paths.promice_data, 'Data', 'artifact_metadata')
   report = icemodel.verification.auditArtifacts( ...
      families="promice", evaluation_data_root=eval_root, ...
      input_data_root=input_root);
   codes = string({report.findings.code});
   gaps = report.findings(codes == "promice_shortwave_darkness_gap");
   testCase.verifyEqual(sort(string({gaps.channel})), ["swd", "swu"]);

   % Location is part of the source metadata contract needed to prove darkness.
   metadata = rmfield(metadata, {'lat', 'lon'});
   Data.Properties.UserData = metadata;
   artifact_metadata = metadata;
   save(paths.promice_data, 'Data', 'artifact_metadata')
   report = icemodel.verification.auditArtifacts( ...
      families="promice", evaluation_data_root=eval_root, ...
      input_data_root=input_root);
   testCase.verifyTrue(any(string({report.findings.code}) ...
      == "promice_shortwave_location_missing"));
   metadata.lat = 70;
   metadata.lon = -40;

   % A non-hourly userdata table is outside the selector's interval contract.
   % The audit must not reinterpret its 15-minute rows as one-hour bins.
   Data.Time = Data.Time(1) + minutes(0:15:45)';
   Data.Properties.UserData = metadata;
   artifact_metadata = metadata;
   save(paths.promice_data, 'Data', 'artifact_metadata')
   report = icemodel.verification.auditArtifacts( ...
      families="promice", evaluation_data_root=eval_root, ...
      input_data_root=input_root);
   testCase.verifyFalse(any(string({report.findings.code}) ...
      == "promice_shortwave_darkness_gap"));

   % One sample cannot establish the hourly interval contract and is skipped.
   Data = Data(1, :);
   save(paths.promice_data, 'Data', 'artifact_metadata')
   report = icemodel.verification.auditArtifacts( ...
      families="promice", evaluation_data_root=eval_root, ...
      input_data_root=input_root);
   testCase.verifyFalse(any(string({report.findings.code}) ...
      == "promice_shortwave_darkness_gap"));

   % Restoring the physical zeros clears the source check without changing the
   % canonical Time row-dimension or station provenance.
   Data = hourly_Data;
   Data.swd = zeros(height(Data), 1);
   Data.swu = zeros(height(Data), 1);
   Data.Properties.UserData = metadata;
   artifact_metadata = metadata;
   save(paths.promice_data, 'Data', 'artifact_metadata')
   report = icemodel.verification.auditArtifacts( ...
      families="promice", evaluation_data_root=eval_root, ...
      input_data_root=input_root);
   testCase.verifyFalse(any(string({report.findings.code}) ...
      == "promice_shortwave_darkness_gap"));
   testCase.verifyEqual(string(Data.Properties.DimensionNames{1}), "Time");

   % A 15-minute met gap at the same deep-night date must not be interpreted as
   % a missing whole-hour zero by this native-userdata policy.
   loaded = load(paths.promice_met, 'met', 'artifact_metadata');
   met = loaded.met;
   met.Properties.DimensionNames{1} = 'Time';
   met.swd(2) = NaN;
   metadata = met.Properties.UserData;
   metadata.met_resample_expected_missing_counts = timetableMissingCounts(met);
   met.Properties.UserData = metadata;
   artifact_metadata = metadata;
   save(paths.promice_met, 'met', 'artifact_metadata')
   report = icemodel.verification.auditArtifacts( ...
      families="promice", evaluation_data_root=eval_root, ...
      input_data_root=input_root);
   testCase.verifyFalse(any(string({report.findings.code}) ...
      == "promice_shortwave_darkness_gap"));
end

function test_mar_reduced_source_fallback_warns_without_daily_constancy(testCase)
   % Reduced-source fixtures may explicitly retain hourly RUH/SMBH. That durable
   % fallback is visible to QA but must not be misdiagnosed as daily-rate drift.
   [eval_root, input_root, paths] = writeAuditTree(testCase.TestData.tmp);
   loaded = load(paths.mar_data, 'Data', 'artifact_metadata');
   loaded.Data.runoff = linspace(0, 0.02, height(loaded.Data))';
   loaded.Data.smb = linspace(-0.01, 0.01, height(loaded.Data))';
   loaded.Data.Properties.UserData.mar_qc_status = 'not_applicable';
   loaded.Data.Properties.UserData.mar_qc_fallback = ...
      'hourly_RUH_SMBH_retained_native_daily_unavailable';
   loaded.Data.Properties.UserData.mar_qc_channels = strings(0, 1);
   loaded.Data.Properties.UserData.mar_qc_replaced_runoff_count = 0;
   loaded.Data.Properties.UserData.mar_qc_replaced_smb_count = 0;
   ndays = numel(loaded.Data.Properties.UserData.mar_qc_runoff_day_status);
   loaded.Data.Properties.UserData.mar_qc_runoff_day_status = ...
      uint8(3 * ones(ndays, 1));
   loaded.Data.Properties.UserData.mar_qc_smb_day_status = ...
      uint8(3 * ones(ndays, 1));
   loaded.Data.Properties.UserData.mar_qc_runoff_daily_reference_mwe = ...
      nan(ndays, 1);
   loaded.Data.Properties.UserData.mar_qc_smb_daily_reference_mwe = ...
      nan(ndays, 1);
   for channel = ["runoff", "smb"]
      loaded.Data.Properties.UserData.( ...
         "mar_qc_preserved_" + channel + "_day_count") = 0;
      loaded.Data.Properties.UserData.( ...
         "mar_qc_replaced_" + channel + "_day_count") = 0;
      loaded.Data.Properties.UserData.( ...
         "mar_qc_unverified_" + channel + "_day_count") = ndays;
   end
   loaded.Data.Properties.UserData.mar_qc_basis = ...
      'MAR native daily RU/SMB unavailable; retained hourly RUH/SMBH';
   Data = loaded.Data;
   artifact_metadata = loaded.Data.Properties.UserData;
   save(paths.mar_data, 'Data', 'artifact_metadata')

   report = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="promice");
   codes = string({report.findings.code});

   testCase.verifyTrue(report.summary.passed);
   testCase.verifyEqual(nnz(codes == "mar_native_daily_unavailable"), 1);
   testCase.verifyFalse(any(codes == "mar_daily_constraint_inconsistent"));
end

function test_mar_ledger_semantics_block_repair_currentness(testCase)
   % A finite reference cannot be called unverified on a complete day, and
   % preserved hourly structure cannot survive a partial boundary. Audit rejects
   % both; bounded alignment can repair only the deterministic boundary defect.
   [eval_root, input_root, paths] = writeAuditTree(testCase.TestData.tmp);
   original = load(paths.mar_data, 'Data', 'artifact_metadata');

   % Contradict one finite daily reference with the unverified status code.
   loaded = original;
   metadata = loaded.Data.Properties.UserData;
   complete_index = 2;
   old_status = metadata.mar_qc_runoff_day_status(complete_index);
   metadata.mar_qc_runoff_day_status(complete_index) = uint8(3);
   if old_status == 1
      metadata.mar_qc_preserved_runoff_day_count = ...
         metadata.mar_qc_preserved_runoff_day_count - 1;
   else
      metadata.mar_qc_replaced_runoff_day_count = ...
         metadata.mar_qc_replaced_runoff_day_count - 1;
   end
   metadata.mar_qc_unverified_runoff_day_count = ...
      metadata.mar_qc_unverified_runoff_day_count + 1;
   loaded.Data.Properties.UserData = metadata;
   Data = loaded.Data;
   artifact_metadata = metadata;
   save(paths.mar_data, 'Data', 'artifact_metadata')

   report = icemodel.verification.auditArtifacts( ...
      families="promice", evaluation_data_root=eval_root, ...
      input_data_root=input_root);
   testCase.verifyTrue(any(string({report.findings.code}) ...
      == "mar_qc_day_ledger_invalid"));
   repair = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice");
   record = repair.records(contains(string({repair.records.filename}), ...
      "site_mar3.11"));
   testCase.verifyEqual(string(record.mar_qc_status), "not_requested");

   % The fixture's first day is partial and therefore cannot be preserved.
   loaded = original;
   metadata = loaded.Data.Properties.UserData;
   testCase.assertEqual(metadata.mar_qc_runoff_day_status(1), uint8(2));
   metadata.mar_qc_runoff_day_status(1) = uint8(1);
   metadata.mar_qc_replaced_runoff_day_count = ...
      metadata.mar_qc_replaced_runoff_day_count - 1;
   metadata.mar_qc_preserved_runoff_day_count = ...
      metadata.mar_qc_preserved_runoff_day_count + 1;
   loaded.Data.Properties.UserData = metadata;
   Data = loaded.Data;
   artifact_metadata = metadata;
   save(paths.mar_data, 'Data', 'artifact_metadata')

   report = icemodel.verification.auditArtifacts( ...
      families="promice", evaluation_data_root=eval_root, ...
      input_data_root=input_root);
   testCase.verifyTrue(any(string({report.findings.code}) ...
      == "mar_qc_day_ledger_invalid"));
   repair = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice");
   record = repair.records(contains(string({repair.records.filename}), ...
      "site_mar3.11"));
   testCase.verifyEqual(string(record.mar_qc_status), "current");
   testCase.verifyTrue(any(string(record.actions) == ...
      "align_mar_daily_metadata"));
end

function test_unverified_merra_orientation_is_an_explicit_blocker(testCase)
   % Numeric MERRA fluxes without a durable orientation marker are unresolved
   % source-quality blockers, not silently accepted warnings.
   [eval_root, input_root, paths] = writeAuditTree(testCase.TestData.tmp);
   loaded = load(paths.merra_data, 'Data', 'artifact_metadata');
   loaded.Data.Properties.UserData = rmfield(loaded.Data.Properties.UserData, ...
      'merra_flux_sign_convention');
   Data = loaded.Data;
   artifact_metadata = loaded.Data.Properties.UserData;
   save(paths.merra_data, 'Data', 'artifact_metadata')

   report = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="promice");
   codes = string({report.findings.code});
   blocker = report.findings(codes == "merra_sign_unverified");

   testCase.verifyFalse(report.summary.passed);
   testCase.verifyEqual(report.summary.blocker_count, 1);
   testCase.verifyEqual(string(blocker.severity), "blocker");
end

function test_merra_time_semantics_are_required_and_axis_checked(testCase)
   % QA distinguishes an invalid policy marker from a leaked native-center axis.
   [eval_root, input_root, paths] = writeAuditTree(testCase.TestData.tmp);
   loaded = load(paths.merra_data, 'Data', 'artifact_metadata');
   loaded.Data.Properties.UserData.merra_time_relabel_policy = 'native_center';
   loaded.Data.Time = loaded.Data.Time + minutes(30);
   Data = loaded.Data;
   artifact_metadata = Data.Properties.UserData;
   save(paths.merra_data, 'Data', 'artifact_metadata')

   report = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="promice");
   codes = string({report.findings.code});

   verifyCodes(testCase, codes, ["merra_time_semantics_invalid", ...
      "merra_time_axis_not_interval_start"]);
end

function test_reduced_merra_artifact_still_requires_time_semantics(testCase)
   % Glacier-only MERRA artifacts have no turbulent-flux sign contract, but their
   % tavg3 support still requires the same durable timestamp provenance.
   [eval_root, input_root, paths] = writeAuditTree(testCase.TestData.tmp);
   loaded = load(paths.merra_data, 'Data', 'artifact_metadata');
   reduced = timetable(loaded.Data.Time, ones(height(loaded.Data), 1), ...
      VariableNames="runoff");
   metadata = rmfield(loaded.artifact_metadata, ...
      ["merra_source_time_coordinate", "merra_time_relabel_policy", ...
      "merra_time_upsample_policy", "merra_collection_support_hours"]);
   reduced.Properties.UserData = metadata;
   Data = reduced;
   artifact_metadata = metadata;
   save(paths.merra_data, 'Data', 'artifact_metadata')

   report = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="promice");
   codes = string({report.findings.code});

   testCase.verifyTrue(any(codes == "merra_time_semantics_missing"));
   testCase.verifyFalse(any(codes == "merra_sign_unverified"));
end

function test_merra_tavg3_ramp_fails_despite_canonical_markers(testCase)
   % First-class QA checks values as well as markers so stale provenance cannot
   % hide an hourly ramp through a three-hour glacier-collection support block.
   [eval_root, input_root, paths] = writeAuditTree(testCase.TestData.tmp);
   loaded = load(paths.merra_data, 'Data', 'artifact_metadata');
   loaded.Data.runoff = (1:height(loaded.Data))';
   Data = loaded.Data;
   artifact_metadata = Data.Properties.UserData;
   save(paths.merra_data, 'Data', 'artifact_metadata')

   report = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="promice");
   codes = string({report.findings.code});

   testCase.verifyTrue(any(codes == "merra_tavg3_support_invalid"));
end

function test_undocumented_required_placeholder_fails(testCase)
   % An all-NaN required channel is valid only when source metadata explicitly
   % describes it as a placeholder rather than an observation.
   [eval_root, input_root, paths] = writeAuditTree(testCase.TestData.tmp);
   loaded = load(paths.promice_met, 'met', 'artifact_metadata');
   loaded.met.Properties.UserData.precip_policy = 'PROMICE precipitation';
   loaded.artifact_metadata = loaded.met.Properties.UserData;
   met = loaded.met;
   artifact_metadata = loaded.artifact_metadata;
   save(paths.promice_met, 'met', 'artifact_metadata')

   report = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="promice");
   codes = string({report.findings.code});

   verifyCodes(testCase, codes, ["undocumented_all_missing", ...
      "promice_precip_policy_invalid"]);
   ppt_findings = report.findings(codes == "undocumented_all_missing" ...
      & string({report.findings.channel}) == "ppt");
   testCase.verifyEqual(string(ppt_findings.severity), "error");
end

function test_partial_required_met_gap_warns_without_failing(testCase)
   % Partial required-channel gaps remain source-faithful but must be visible as
   % non-runnable forcing rather than buried only in the channel ledger.
   [eval_root, input_root, paths] = writeAuditTree(testCase.TestData.tmp);
   loaded = load(paths.promice_met, 'met');
   loaded.met.albedo(4:5) = NaN;
   loaded.met.Properties.UserData.met_resample_expected_missing_counts.albedo = 2;
   saveMet(paths.promice_met, loaded.met)

   report = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="promice", report_dir="");
   codes = string({report.findings.code});
   findings = report.findings(codes == "required_met_gap" ...
      & string({report.findings.channel}) == "albedo");

   testCase.verifyTrue(report.summary.passed);
   testCase.verifyEqual(numel(findings), 1);
   testCase.verifyEqual(string(findings.severity), "warning");
   testCase.verifyTrue(contains(string(findings.message), "2 of 13"));
   testCase.verifyTrue(contains(string(findings.message), "% finite"));
   testCase.verifyTrue(contains(string(findings.message), ...
      "not directly runnable"));
end

function test_detects_manifest_schema_and_source_list_drift(testCase)
   % Required fields and staged source/file declarations should be checked before
   % valid leftover files can make a stale manifest appear usable.
   [eval_root, input_root, paths] = writeAuditTree(testCase.TestData.tmp);
   manifest = jsondecode(fileread(paths.manifest));
   manifest.cases = rmfield(manifest.cases, 'native_timestep');
   manifest.cases.forcing_sources = {'promice', 'missing_forcing'};
   manifest.cases.eval_sources = [manifest.cases.eval_sources; {'missing_eval'}];
   manifest.cases.colocation.empty_leg = struct('staged', true, ...
      'source_id', 'missing_eval');
   writeJson(paths.manifest, manifest)

   report = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="promice");
   codes = string({report.findings.code});

   verifyCodes(testCase, codes, ["missing_case_field", ...
      "forcing_source_without_met", "eval_source_without_artifact", ...
      "staged_leg_without_artifact"]);
end

function test_atomic_esm_runtime_met_is_audited_exactly(testCase)
   % Atomic ESM manifests omit met paths, so QA must use the same exact runtime
   % selection as plotting and ignore unrelated wildcard-compatible siblings.
   [eval_root, input_root, paths] = writeEsmAuditTree(testCase.TestData.tmp);
   manifest_before = fileBytes(paths.manifest);
   before = fileSnapshot([eval_root; input_root]);

   report = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="esm_snowmip");

   % The nested exact file wins over a decoy and contributes every met channel.
   after = fileSnapshot([eval_root; input_root]);
   testCase.verifyEqual(after, before);
   met_records = report.artifacts(string({report.artifacts.kind}) == "met");
   testCase.verifyNumElements(met_records, 1);
   testCase.verifyEqual(string(met_records.path), string(paths.met_nested));
   testCase.verifyEqual(string(met_records.source), "esm_snowmip");
   testCase.verifyTrue(met_records.exists);
   testCase.verifyEqual(met_records.n_channels, 11);
   testCase.verifyEqual(report.summary.case_count, 1);
   testCase.verifyEqual(report.summary.artifact_count, 3);
   testCase.verifyEqual(report.summary.channel_count, 12);
   testCase.verifyTrue(report.summary.passed);
   testCase.verifyFalse(any(string({report.artifacts.path}) == paths.decoy));

   % The canonical runtime resolver supports the legacy flat layout after the
   % nested exact artifact moves, without admitting the nested decoy.
   movefile(paths.met_nested, paths.met_flat)
   before = fileSnapshot([eval_root; input_root]);
   report = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="esm_snowmip");
   after = fileSnapshot([eval_root; input_root]);
   testCase.verifyEqual(after, before);
   met_records = report.artifacts(string({report.artifacts.kind}) == "met");
   testCase.verifyNumElements(met_records, 1);
   testCase.verifyEqual(string(met_records.path), string(paths.met_flat));
   testCase.verifyEqual(met_records.n_channels, 11);

   % A missing exact artifact remains in the ledger at its canonical flat path
   % and becomes an actionable missing_artifact instead of disappearing.
   delete(paths.met_flat)
   report = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="esm_snowmip");
   met_records = report.artifacts(string({report.artifacts.kind}) == "met");
   testCase.verifyNumElements(met_records, 1);
   testCase.verifyEqual(string(met_records.path), string(paths.met_flat));
   testCase.verifyFalse(met_records.exists);
   testCase.verifyTrue(any(string({report.findings.code}) == "missing_artifact"));

   % A genuine resolver failure is reported rather than silently suppressing
   % the entire forcing leg. The atomic manifest itself remains byte-stable.
   rmdir(input_root, 's')
   report = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="esm_snowmip");
   testCase.verifyTrue(any(string({report.findings.code}) ...
      == "runtime_met_resolution_error"));
   testCase.verifyEqual(fileBytes(paths.manifest), manifest_before);
   manifest = jsondecode(fileread(paths.manifest));
   testCase.verifyFalse(isfield(manifest.cases, 'forcing_sources'));
   testCase.verifyFalse(isfield(manifest.cases, 'colocation'));
   testCase.verifyFalse(isfield(manifest.cases, 'met_files'));
end

function test_repairRcmArtifactMetadata_repairs_racmo_subl_once(testCase)
   % Legacy RACMO sublimation is negative for surface loss. The bounded repair
   % must flip only that channel, stamp both conventions, and be byte-stable on
   % pass two while the artifact audit enforces the durable contract.
   [eval_root, input_root, paths] = writeAuditTree(testCase.TestData.tmp);
   racmo_file = attachRacmoSublFixture(input_root, paths.manifest);
   original = load(racmo_file, 'Data', 'artifact_metadata', 'auxiliary');
   before = fileBytes(racmo_file);

   audit = icemodel.verification.auditArtifacts( ...
      families="promice", evaluation_data_root=eval_root, ...
      input_data_root=input_root);
   testCase.verifyTrue(any(string({audit.findings.code}) ...
      == "racmo_subl_sign_unverified"));

   sign_disabled = ...
      icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice", ...
      source_id="racmo2.3p3", repair_racmo_subl=false);
   testCase.verifyEqual(string(sign_disabled.records.status), "unchanged");
   testCase.verifyFalse(any(contains(string(sign_disabled.records.actions), ...
      "racmo_subl")));

   % The source filter must exclude every manifest reference when no source
   % matches the requested identifier.
   excluded = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice", ...
      source_id="racmo9.9");
   testCase.verifyEmpty(excluded.records);

   dry = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice", ...
      source_id="racmo2.3p3");
   testCase.verifyEqual(string(dry.records.status), "would_repair");
   testCase.verifyEqual(fileBytes(racmo_file), before);

   first = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice", ...
      source_id="racmo2.3p3", dry_run=false);
   testCase.verifyEqual(string(first.records.status), "repaired");
   testCase.verifyTrue(any(string(first.records.actions) ...
      == "flip_racmo_subl_sign"));
   testCase.verifyEqual(string(first.records.changed_variables), "subl");
   expected_metadata = ["racmo_subl_native_sign_convention"; ...
      "racmo_subl_sign_convention"];
   testCase.verifyEqual(sort(string( ...
      first.records.changed_metadata_fields(:))), sort(expected_metadata));
   testCase.verifyTrue(first.records.unrelated_payload_preserved);
   testCase.verifyEqual(strlength(first.records.hash_before), 64);
   testCase.verifyEqual(strlength(first.records.hash_after), 64);
   testCase.verifyNotEqual(first.records.hash_before, first.records.hash_after);
   repaired = load(racmo_file, 'Data', 'artifact_metadata', 'auxiliary');
   testCase.verifyEqual(repaired.Data.Time, original.Data.Time);
   testCase.verifyEqual(repaired.Data.guard, original.Data.guard);
   testCase.verifyEqual(repaired.Data.subl, -original.Data.subl);
   testCase.verifyEqual(repaired.auxiliary, original.auxiliary);
   testCase.verifyEqual(string( ...
      repaired.artifact_metadata.racmo_subl_native_sign_convention), ...
      "negative_loss_positive_deposition");
   testCase.verifyEqual(string( ...
      repaired.Data.Properties.UserData.racmo_subl_sign_convention), ...
      "positive_loss_negative_deposition");

   after_first = fileBytes(racmo_file);
   second = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice", ...
      source_id="racmo2.3p3", dry_run=false);
   testCase.verifyEqual(string(second.records.status), "unchanged");
   testCase.verifyEqual(second.records.hash_before, first.records.hash_after);
   testCase.verifyEqual(second.records.hash_after, first.records.hash_after);
   testCase.verifyEqual(fileBytes(racmo_file), after_first);

   audit = icemodel.verification.auditArtifacts( ...
      families="promice", evaluation_data_root=eval_root, ...
      input_data_root=input_root);
   testCase.verifyFalse(any(startsWith(string({audit.findings.code}), ...
      "racmo_subl_sign")));
end

function test_racmo_subl_repair_rejects_ambiguous_markers(testCase)
   % A partial, unknown, or nonscalar marker cannot prove whether the numeric
   % sign was already changed. The repair must stop instead of double-flipping;
   % RACMO artifacts without subl remain outside this contract.
   [eval_root, input_root, paths] = writeAuditTree(testCase.TestData.tmp);
   racmo_file = attachRacmoSublFixture(input_root, paths.manifest);
   original = load(racmo_file, 'Data', 'artifact_metadata', 'auxiliary');

   loaded = original;
   loaded.artifact_metadata.racmo_subl_native_sign_convention = ...
      'negative_loss_positive_deposition';
   loaded.Data.Properties.UserData = loaded.artifact_metadata;
   Data = loaded.Data;
   artifact_metadata = loaded.artifact_metadata;
   auxiliary = loaded.auxiliary;
   save(racmo_file, 'Data', 'artifact_metadata', 'auxiliary')
   report = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice", ...
      source_id="racmo2.3p3");
   testCase.verifyEqual(string(report.records.status), "restage_required");
   testCase.verifyTrue(contains(string(report.records.reason), "incomplete"));

   loaded.artifact_metadata.racmo_subl_native_sign_convention = ...
      ["negative_loss_positive_deposition", "unknown"];
   loaded.artifact_metadata.racmo_subl_sign_convention = ...
      'positive_loss_negative_deposition';
   loaded.Data.Properties.UserData = loaded.artifact_metadata;
   Data = loaded.Data;
   artifact_metadata = loaded.artifact_metadata;
   save(racmo_file, 'Data', 'artifact_metadata', 'auxiliary')
   report = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice", ...
      source_id="racmo2.3p3");
   testCase.verifyEqual(string(report.records.status), "restage_required");
   testCase.verifyTrue(contains(string(report.records.reason), "scalar"));

   % Recognized-but-wrong scalar metadata is an explicit QA error.
   loaded.artifact_metadata.racmo_subl_native_sign_convention = 'unknown';
   loaded.artifact_metadata.racmo_subl_sign_convention = 'unknown';
   loaded.Data.Properties.UserData = loaded.artifact_metadata;
   Data = loaded.Data;
   artifact_metadata = loaded.artifact_metadata;
   save(racmo_file, 'Data', 'artifact_metadata', 'auxiliary')
   audit = icemodel.verification.auditArtifacts( ...
      families="promice", evaluation_data_root=eval_root, ...
      input_data_root=input_root);
   testCase.verifyTrue(any(string({audit.findings.code}) ...
      == "racmo_subl_sign_not_canonical"));

   % A RACMO file with no subl channel is ignored by both repair and QA.
   Data = removevars(original.Data, "subl");
   artifact_metadata = original.artifact_metadata;
   Data.Properties.UserData = artifact_metadata;
   save(racmo_file, 'Data', 'artifact_metadata', 'auxiliary')
   report = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice", ...
      source_id="racmo2.3p3", dry_run=false);
   testCase.verifyTrue(ismember(string(report.records.status), ...
      ["repaired", "unchanged"]));
   testCase.verifyFalse(any(contains(string(report.records.actions), ...
      "racmo_subl")));
   second = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice", ...
      source_id="racmo2.3p3", dry_run=false);
   testCase.verifyEqual(string(second.records.status), "unchanged");
   audit = icemodel.verification.auditArtifacts( ...
      families="promice", evaluation_data_root=eval_root, ...
      input_data_root=input_root);
   testCase.verifyFalse(any(startsWith(string({audit.findings.code}), ...
      "racmo_subl_sign")));
end

function test_repairMetTimeSupport_dry_write_and_idempotence(testCase)
   % Legacy source rows are recoverable without raw data; dry-run is byte-stable.
   [~, ~, paths] = writeAuditTree(testCase.TestData.tmp);
   loaded = load(paths.promice_met, 'met');
   source = loaded.met(1:4:end, :);
   source.swd = (0:100:300)';
   legacy_time = (source.Time(1):minutes(15):source.Time(end))';
   met = retime(source, legacy_time, 'linear');
   metadata = source.Properties.UserData;
   metadata.met_resample_policy = "linear_adjacent_finite_only";
   metadata.met_resample_source_row_count = height(source);
   metadata.met_resample_source_cadence_seconds = 3600;
   metadata.met_resample_source_time_gap_count = 0;
   metadata.met_resample_source_missing_counts = timetableMissingCounts(source);
   metadata.met_resample_expected_missing_counts = timetableMissingCounts(met);
   met.Properties.UserData = metadata;
   auxiliary = struct('label', "preserve", 'values', [1, 2, 3]);
   save(paths.promice_met, 'met', 'auxiliary')
   before = fileBytes(paths.promice_met);

   dry = icemodel.verification.setup.repairMetTimeSupport( ...
      paths.promice_met);
   testCase.verifyEqual(dry.summary.would_repair_count, 1);
   testCase.verifyEqual(fileBytes(paths.promice_met), before);

   written = icemodel.verification.setup.repairMetTimeSupport( ...
      paths.promice_met, dry_run=false);
   repaired = load(paths.promice_met, 'met', 'artifact_metadata', 'auxiliary');
   testCase.verifyEqual(written.summary.repaired_count, 1);
   testCase.verifyEqual(height(repaired.met), 4 * height(source));
   testCase.verifyEqual(repaired.met.Time(end), ...
      source.Time(end) + minutes(45));
   testCase.verifyEqual(repaired.met.swd, repelem(source.swd, 4));
   testCase.verifyEqual(string( ...
      repaired.artifact_metadata.met_resample_policy), ...
      "interval_start_zero_order_hold");
   testCase.verifyEqual(repaired.auxiliary, auxiliary);

   second = icemodel.verification.setup.repairMetTimeSupport( ...
      paths.promice_met);
   testCase.verifyEqual(second.summary.unchanged_count, 1);
   testCase.verifyEqual(second.summary.would_repair_count, 0);
end

function test_repairMetTimeSupport_rejects_missing_and_unknown(testCase)
   % Exact repair inputs fail closed when a path or provenance policy is unknown.
   [~, ~, paths] = writeAuditTree(testCase.TestData.tmp);
   testCase.verifyError(@() ...
      icemodel.verification.setup.repairMetTimeSupport( ...
      fullfile(testCase.TestData.tmp, 'missing.mat')), ...
      'icemodel:verification:repairMetTimeSupport:fileNotFound');

   loaded = load(paths.promice_met, 'met', 'artifact_metadata');
   loaded.met.Properties.UserData.met_resample_policy = "unknown";
   met = loaded.met;
   artifact_metadata = met.Properties.UserData;
   save(paths.promice_met, 'met', 'artifact_metadata')
   testCase.verifyError(@() ...
      icemodel.verification.setup.repairMetTimeSupport(paths.promice_met), ...
      'icemodel:verification:repairMetTimeSupport:unsupportedPolicy');
end

function test_repairMetTimeSupport_rejects_invented_omitted_source_row(testCase)
   % A legacy regular 15-minute grid can contain a finite row at an omitted
   % native timestamp. More recovered grid rows than recorded source rows is not
   % sufficient evidence for an exact source-light repair.
   [~, ~, paths] = writeAuditTree(testCase.TestData.tmp);
   loaded = load(paths.promice_met, 'met');
   met = loaded.met(1:9, :);
   met.swd = (0:8)';
   metadata = met.Properties.UserData;
   metadata.met_resample_policy = "linear_adjacent_finite_only";
   metadata.met_resample_source_row_count = 2;
   metadata.met_resample_source_cadence_seconds = 3600;
   metadata.met_resample_source_time_gap_count = 0;
   met.Properties.UserData = metadata;
   artifact_metadata = metadata;
   save(paths.promice_met, 'met', 'artifact_metadata')

   testCase.verifyError(@() ...
      icemodel.verification.setup.repairMetTimeSupport(paths.promice_met), ...
      'icemodel:verification:repairMetTimeSupport:ambiguousLegacyGrid');
end

function test_repairMetTimeSupport_rejects_recorded_source_gap(testCase)
   % Even a row-count-compatible legacy payload is unrecoverable when provenance
   % says the native source omitted a cadence interval.
   [~, ~, paths] = writeAuditTree(testCase.TestData.tmp);
   loaded = load(paths.promice_met, 'met');
   source = loaded.met(1:4:end, :);
   met = retime(source, (source.Time(1):minutes(15):source.Time(end))', ...
      'linear');
   metadata = source.Properties.UserData;
   metadata.met_resample_policy = "linear_adjacent_finite_only";
   metadata.met_resample_source_row_count = height(source);
   metadata.met_resample_source_cadence_seconds = 3600;
   metadata.met_resample_source_time_gap_count = 1;
   met.Properties.UserData = metadata;
   artifact_metadata = metadata;
   save(paths.promice_met, 'met', 'artifact_metadata')

   testCase.verifyError(@() ...
      icemodel.verification.setup.repairMetTimeSupport(paths.promice_met), ...
      'icemodel:verification:repairMetTimeSupport:ambiguousLegacyGrid');
end

function test_repairMetTimeSupport_rejects_unproven_merra_tavg3(testCase)
   % Recovering hourly rows from legacy MERRA met is unsafe when no native glc
   % inventory proves that the 00/03/... rows were not themselves interpolated.
   filename = fullfile(testCase.TestData.tmp, ...
      'met_site_merra2_20000101_20000101_15m.mat');
   writeLegacyMerraMet(filename, false);

   testCase.verifyError(@() ...
      icemodel.verification.setup.repairMetTimeSupport(filename), ...
      'icemodel:verification:repairMetTimeSupport:unprovenMerraSourceGrid');
end

function test_repairMetTimeSupport_detects_legacy_merra_label(testCase)
   % The old `_merra_` storage token must route through MERRA safety even when
   % its legacy UserData contains only generic met-resampling provenance.
   filename = fullfile(testCase.TestData.tmp, ...
      'met_site_merra_20000101_20000101_15m.mat');
   writeLegacyMerraMet(filename, false, false);

   testCase.verifyError(@() ...
      icemodel.verification.setup.repairMetTimeSupport(filename), ...
      'icemodel:verification:repairMetTimeSupport:unprovenMerraSourceGrid');
end

function test_repairMetTimeSupport_reholds_proven_merra_tavg3(testCase)
   % With an exact raw-grid proof, repair first restores each three-hour glc hold
   % on hourly rows and only then repeats those hours onto the 15-minute grid.
   filename = fullfile(testCase.TestData.tmp, ...
      'met_site_merra2_20000101_20000101_15m.mat');
   writeLegacyMerraMet(filename, true);

   report = icemodel.verification.setup.repairMetTimeSupport( ...
      filename, dry_run=false);
   repaired = load(filename, 'met', 'artifact_metadata');

   testCase.verifyEqual(report.summary.repaired_count, 1);
   testCase.verifyEqual(repaired.met.runoff, ...
      [ones(12, 1); 4 * ones(4, 1)]);
   testCase.verifyTrue( ...
      icemodel.forcing.helpers.hasConstantMerraTavg3Support(repaired.met));
   testCase.verifyEqual(repaired.met.Properties.UserData, ...
      repaired.artifact_metadata);
end

function test_repairRcmArtifactMetadata_accepts_current_15m_merra_met(testCase)
   % Current 15-minute met has already crossed the resampling boundary; metadata
   % repair must validate it without sending it through the hourly tavg3 helper.
   [eval_root, input_root, paths] = writeAuditTree(testCase.TestData.tmp);
   met_file = attachMerraMetFixture(input_root, paths.manifest, ...
      "interval_start_zero_order_hold");
   before = fileBytes(met_file);

   report = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice");
   filenames = string({report.records.filename});
   is_met = endsWith(filenames, ...
      "met_site_merra2_20000101_20000101_15m.mat");
   testCase.assertEqual(nnz(is_met), 1, strjoin(filenames, newline));
   record = report.records(is_met);

   testCase.verifyEqual(string(record.status), "unchanged");
   testCase.verifyEmpty(record.actions);
   testCase.verifyEqual(fileBytes(met_file), before);
end

function test_repairRcmArtifactMetadata_proves_native_merra_grid_once(testCase)
   % A legacy hourly ramp is repairable only after the native daily glc time axis
   % proves every saved 3-hour source row; the durable proof makes pass two light.
   [eval_root, input_root, paths] = writeAuditTree(testCase.TestData.tmp);
   merra_dir = fullfile(testCase.TestData.tmp, 'merra-source');
   writeTinyMerraGlcSource(merra_dir, datetime(2000, 1, 1, TimeZone="UTC"));
   loaded = load(paths.merra_data, 'Data', 'artifact_metadata');
   loaded.Data.runoff = (1:height(loaded.Data))';
   loaded.Data.Properties.UserData = loaded.artifact_metadata;
   Data = loaded.Data;
   artifact_metadata = loaded.artifact_metadata;
   auxiliary = struct('label', "preserve", 'values', [1, 2, 3]);
   save(paths.merra_data, 'Data', 'artifact_metadata', 'auxiliary')

   first = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice", ...
      merra_dir=string(merra_dir), dry_run=false);
   repaired = load(paths.merra_data, ...
      'Data', 'artifact_metadata', 'auxiliary');
   record = first.records(string({first.records.source_id}) == "merra2");

   testCase.verifyEqual(string(record.status), "repaired");
   testCase.verifyTrue(any(string(record.actions) ...
      == "stamp_merra_tavg3_source_grid"));
   testCase.verifyEqual(first.merra_source_reads, 1);
   testCase.verifyEqual(repaired.Data.runoff, [1; 1; 1; 4]);
   testCase.verifyEqual(repaired.auxiliary, auxiliary);
   testCase.verifyTrue( ...
      icemodel.forcing.helpers.hasProvenMerraTavg3SourceGrid( ...
      repaired.Data, repaired.artifact_metadata));

   second = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice", ...
      merra_dir=string(merra_dir));
   record = second.records(string({second.records.source_id}) == "merra2");
   testCase.verifyEqual(string(record.status), "unchanged");
   testCase.verifyEqual(second.merra_source_reads, 0);
end

function test_repairRcmArtifactMetadata_rejects_legacy_15m_merra_met(testCase)
   % A linear-era 15-minute MERRA artifact cannot be treated as hourly Data; it
   % must remain byte-stable and route to the dedicated met repair/restage path.
   [eval_root, input_root, paths] = writeAuditTree(testCase.TestData.tmp);
   met_file = attachMerraMetFixture(input_root, paths.manifest, ...
      "linear_adjacent_finite_only");
   before = fileBytes(met_file);

   report = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice");
   filenames = string({report.records.filename});
   is_met = endsWith(filenames, ...
      "met_site_merra2_20000101_20000101_15m.mat");
   testCase.assertEqual(nnz(is_met), 1, strjoin(filenames, newline));
   record = report.records(is_met);

   testCase.verifyEqual(string(record.status), "restage_required");
   testCase.verifyTrue(contains(string(record.reason), "repairMetTimeSupport"));
   testCase.verifyEqual(fileBytes(met_file), before);
end

function test_repairRcmArtifactMetadata_rejects_bridged_15m_merra_met(testCase)
   % Canonical-looking markers cannot override the source-derived missing-value
   % lower bound when a 15-minute artifact has bridged omitted support.
   [eval_root, input_root, paths] = writeAuditTree(testCase.TestData.tmp);
   met_file = attachMerraMetFixture(input_root, paths.manifest, ...
      "interval_start_zero_order_hold");
   loaded = load(met_file, 'met', 'artifact_metadata');
   loaded.artifact_metadata.met_resample_expected_missing_counts.runoff = 1;
   loaded.met.Properties.UserData = loaded.artifact_metadata;
   met = loaded.met;
   artifact_metadata = loaded.artifact_metadata;
   save(met_file, 'met', 'artifact_metadata')

   report = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice");
   filenames = string({report.records.filename});
   record = report.records(endsWith(filenames, ...
      "met_site_merra2_20000101_20000101_15m.mat"));

   testCase.verifyEqual(string(record.status), "restage_required");
   testCase.verifyTrue(contains(string(record.reason), "bridges"));
end

function [eval_root, input_root, protocol_file] = writeRetmipAuditTree(root)
   %WRITERETMIPAUDITTREE Build one protocol NetCDF plus status-only native leg.
   eval_root = fullfile(root, "retmip_eval");
   input_root = fullfile(root, "retmip_input");
   family_root = fullfile(eval_root, "retmip");
   case_root = fullfile(family_root, "fixture");
   for folder = [input_root, family_root, case_root]
      mkdir(folder)
   end

   % The top-level evaluation target remains a normal MAT bundle.
   times = (datetime(2000, 1, 1, TimeZone="UTC"):hours(3): ...
      datetime(2000, 1, 1, 6, 0, 0, TimeZone="UTC"))';
   obs = timetable(times, 260 + zeros(numel(times), 1), ...
      VariableNames="tsfc");
   obs = icemodel.forcing.helpers.stampMetadata(obs);
   targets = struct('format', 'timeseries', 'data', obs);
   observations_file = fullfile(case_root, "observations.mat");
   save(observations_file, 'targets')

   % Protocol model output is deliberately native NetCDF.
   protocol_file = string(fullfile(case_root, ...
      "RetMIP_fixture_columns.nc"));
   writeRetmipProtocolNetcdf(protocol_file, "valid");
   colocation = struct();
   colocation.retmip = struct( ...
      'kind', 'protocol_userdata_and_native_met', ...
      'staged', true, ...
      'model_output_files', char(protocol_file), ...
      'model_output_variables', {{'rho'; 'time'}});
   colocation.native_met = struct( ...
      'staged', true, 'status', 'staged', ...
      'source_family', 'fixture', 'source_id', 'SITE');

   case_entry = struct( ...
      'case_id', 'fixture', ...
      'case_type', 'firn_observational', ...
      'site_id', 'SITE', ...
      'site_name', 'RetMIP fixture', ...
      'surface_zone', 'accumulation', ...
      'eval_target', 'surface_temperature', ...
      'permafrost_zone', 'none', ...
      'site_location', struct('lat', 70, 'lon', -40, 'elev', 2000, ...
      'lat_wgs84', 70, 'lon_wgs84', -40), ...
      'period', struct('start', '2000-01-01 00:00:00', ...
      'end', '2000-01-02 00:00:00'), ...
      'evaluation_file', 'fixture/observations.mat', ...
      'forcing_sources', {cell(0, 1)}, ...
      'eval_sources', {{'retmip_protocol'}}, ...
      'comparison_variables', {{'tsfc'}}, ...
      'observation_variables', struct(), ...
      'colocation', colocation, ...
      'native_timestep', '3hr', ...
      'notes', 'RetMIP audit fixture');
   manifest = struct( ...
      'dataset_family', 'retmip', ...
      'source_doi', '', ...
      'source_url', 'https://example.invalid', ...
      'source_version', 'fixture', ...
      'retrieval_date', '2000-01-02', ...
      'cases', case_entry);
   writeJson(fullfile(family_root, "manifest.json"), manifest)
end

function [data_file, flagged] = attachAlbedoTransientFixture( ...
      eval_root, input_root)
   %ATTACHALBEDOTRANSIENTFIXTURE Add one masked IMAU userdata source leg.
   days_utc = (datetime(2017, 1, 1, TimeZone="UTC") ...
      + caldays(0:89))';
   Time = (days_utc(1):hours(1):(days_utc(end) + hours(23)))';
   day = dateshift(Time, 'start', 'day');
   swd = 100 * ones(size(Time));
   alpha = 0.85 * ones(size(Time));
   alpha(day == days_utc(44)) = 0.72;
   alpha(day == days_utc(46)) = 0.40;
   alpha(day == days_utc(48)) = 0.72;
   swu = alpha .* swd;
   swn = swd - swu;
   netr = swn - 50;

   % Preserve the raw shortwave pair while masking every derived channel on
   % the shared detector's recovered episode.
   [flagged, qc] = ...
      icemodel.forcing.helpers.dailyAlbedoAnomalyFlags(Time, swd, swu);
   qc = rmfield(qc, 'diagnostics');
   alpha(flagged) = NaN;
   swn(flagged) = NaN;
   netr(flagged) = NaN;
   Data = timetable(Time, swd, swu, alpha, swn, netr, ...
      VariableNames={'swd', 'swu', 'albedo', 'swn', 'netr'});
   Data = icemodel.forcing.helpers.stampMetadata(Data);
   Data.Properties.UserData = struct( ...
      'source_family', 'imau', ...
      'station', 'S21', ...
      'albedo_transient_qc', qc);

   userdata_root = fullfile(input_root, 'userdata', 'imau');
   mkdir(userdata_root)
   filename = "fixture_imau_20170101_20170331.mat";
   data_file = string(fullfile(userdata_root, filename));
   saveData(data_file, Data)

   % Register the new source through the same colocation contract used by
   % production RetMIP manifests.
   manifest_file = fullfile(eval_root, 'retmip', 'manifest.json');
   manifest = jsondecode(fileread(manifest_file));
   manifest.cases.period.end = '2017-03-31 23:00:00';
   manifest.cases.eval_sources = {'retmip_protocol'; 'imau'};
   manifest.cases.colocation.imau = struct( ...
      'kind', 'station_eval', ...
      'staged', true, ...
      'source', 'imau', ...
      'source_id', 'imau', ...
      'data_files', char(fullfile('imau', filename)), ...
      'window', struct( ...
      'start', '2017-01-01 00:00:00', ...
      'end', '2017-03-31 23:00:00'));
   writeJson(manifest_file, manifest)
end

function writeRetmipProtocolNetcdf(filename, mode)
   %WRITERETMIPPROTOCOLNETCDF Write focused native-format contract variants.
   if isfile(filename)
      delete(filename)
   end
   if mode == "missing_time"
      nccreate(filename, 'rho', Dimensions={'sample', 3});
      ncwrite(filename, 'rho', [300; 350; 400]);
      return
   elseif mode == "missing_variables"
      nccreate(filename, 'unknown', Dimensions={'sample', 3});
      ncwrite(filename, 'unknown', [1; 2; 3]);
      return
   end

   nccreate(filename, 'time', Dimensions={'time', 3});
   nccreate(filename, 'rho', Dimensions={'time', 3});
   ncwrite(filename, 'rho', [300; 350; 400]);
   if mode ~= "missing_units"
      ncwriteatt(filename, 'time', 'Units', ...
         'Days since 2000-01-01 00:00');
   end
   if mode == "nonmonotonic"
      time = [0; 0; 0.25];
   elseif mode == "nonfinite"
      time = [0; NaN; 0.25];
   else
      time = [0; 0.125; 0.25];
   end
   ncwrite(filename, 'time', time);
end

function [eval_root, input_root, observations_file] = writeProfileAuditTree(root)
   %WRITEPROFILEAUDITTREE Build one SUMup-style nested event/profile bundle.
   eval_root = fullfile(root, "profile_eval");
   input_root = fullfile(root, "profile_input");
   family_root = fullfile(eval_root, "sumup");
   case_root = fullfile(family_root, "fixture");
   for folder = [input_root, family_root, case_root]
      mkdir(folder)
   end

   % Survey records intentionally repeat dates across depth rows and are not
   % sorted globally; uncertainty units match their measurement channels.
   survey_time = [ ...
      datetime(2001, 1, 1, TimeZone="UTC")
      datetime(2000, 1, 1, TimeZone="UTC")
      datetime(2000, 1, 1, TimeZone="UTC")];
   density = table([10; 20; 20], [300; 400; 450], [0; 0; 1], [1; 1; 2], ...
      ["d1"; "d2"; "d2"], survey_time, VariableNames={ ...
      'error', 'density', 'start_depth', 'stop_depth', ...
      'measurement_id', 'datetime'});
   density.Properties.VariableUnits = ...
      {'kg m-3', 'kg m-3', 'm', 'm', '', ''};

   temperature = timetable([1; 1; 1], [-20; -19; -18], [0; 0; 1], ...
      ["t1"; "t2"; "t2"], RowTimes=survey_time, VariableNames={ ...
      'error', 'subsurface_temperature', 'depth', 'measurement_id'});
   temperature.Properties.VariableUnits = {'degC', 'degC', 'm', ''};

   start_date = NaT(2, 1, TimeZone="UTC");
   end_date = NaT(2, 1, TimeZone="UTC");
   start_date(2) = datetime(2000, 1, 1, TimeZone="UTC");
   end_date(2) = datetime(2001, 1, 1, TimeZone="UTC");
   smb = table([0.01; 0.02], [0.2; 0.3], start_date, end_date, ...
      [1999; 2000], [2000; 2001], ["s1"; "s2"], VariableNames={ ...
      'error', 'smb', 'start_date', 'end_date', 'start_year', 'end_year', ...
      'measurement_id'});
   smb.Properties.VariableUnits = ...
      {'m w.e.', 'm w.e.', '', '', '1', '1', ''};

   data = struct('format', 'subsurface_profile_bundle', ...
      'density', density, 'subsurface_temperature', temperature, 'smb', smb);
   targets = struct('format', 'subsurface_profile_bundle', ...
      'data', data, 'metadata', struct());
   observations_file = string(fullfile(case_root, "observations.mat"));
   save(observations_file, 'targets')

   colocation = struct('sumup', struct( ...
      'kind', 'subsurface_profile', 'staged', true, ...
      'obs_file', 'fixture/observations.mat'));
   case_entry = struct( ...
      'case_id', 'fixture', ...
      'case_type', 'firn_observational', ...
      'site_id', 'SITE', ...
      'site_name', 'SUMup fixture', ...
      'surface_zone', 'accumulation', ...
      'eval_target', 'subsurface_profile_bundle', ...
      'permafrost_zone', 'none', ...
      'site_location', struct('lat', 70, 'lon', -40, 'elev', 2000, ...
      'lat_wgs84', 70, 'lon_wgs84', -40), ...
      'period', struct('start', '1999-01-01 00:00:00', ...
      'end', '2002-01-01 00:00:00'), ...
      'evaluation_file', 'fixture/observations.mat', ...
      'forcing_sources', {cell(0, 1)}, ...
      'eval_sources', {{'sumup_obs'}}, ...
      'comparison_variables', {{'density'; 'subsurface_temperature'; 'smb'}}, ...
      'observation_variables', struct(), ...
      'colocation', colocation, ...
      'native_timestep', 'event', ...
      'notes', 'SUMup audit fixture');
   manifest = struct( ...
      'dataset_family', 'sumup', ...
      'source_doi', '', ...
      'source_url', 'https://example.invalid', ...
      'source_version', 'fixture', ...
      'retrieval_date', '2002-01-02', ...
      'cases', case_entry);
   writeJson(fullfile(family_root, "manifest.json"), manifest)
end

function [eval_root, input_root, paths] = writeAuditTree(root)
   %WRITEAUDITTREE Build one valid compact PROMICE/MAR/MERRA artifact graph.
   eval_root = fullfile(root, "eval");
   input_root = fullfile(root, "input");
   family_root = fullfile(eval_root, "promice");
   case_root = fullfile(family_root, "site");
   met_root = fullfile(input_root, "met", "promice");
   promice_root = fullfile(input_root, "userdata", "promice");
   mar_root = fullfile(input_root, "userdata", "mar3.11");
   merra_root = fullfile(input_root, "userdata", "merra2");
   for folder = [case_root, met_root, promice_root, mar_root, merra_root]
      mkdir(folder)
   end

   % Case-level observation target: one canonical temperature timeseries.
   short_time = (datetime(2000, 1, 1, 0, 0, 0, TimeZone="UTC"):hours(1): ...
      datetime(2000, 1, 1, 3, 0, 0, TimeZone="UTC"))';
   obs = timetable(short_time, 260 + zeros(numel(short_time), 1), ...
      VariableNames="tsfc");
   obs = icemodel.forcing.helpers.stampMetadata(obs);
   targets = struct('format', 'timeseries', 'data', obs);
   paths.observations = fullfile(case_root, "observations.mat");
   save(paths.observations, 'targets')

   % PROMICE native Data remains hourly and does not invent precipitation.
   Data = timetable(short_time, 260 + zeros(numel(short_time), 1), ...
      VariableNames="tair");
   Data = icemodel.forcing.helpers.stampMetadata(Data);
   Data.Properties.UserData = struct('source', 'promice');
   paths.promice_data = fullfile(promice_root, ...
      "site_promice_20000101_20000101.mat");
   saveData(paths.promice_data, Data)

   % PROMICE met uses the canonical required contract; ppt/snowf are explicit
   % all-NaN placeholders while rainf is source-faithful liquid precipitation.
   met_time = (short_time(1):minutes(15):short_time(end))';
   n = numel(met_time);
   met = timetable(met_time, 260 + zeros(n, 1), 100 + zeros(n, 1), ...
      200 + zeros(n, 1), 0.7 + zeros(n, 1), 5 + zeros(n, 1), ...
      80 + zeros(n, 1), 80000 + zeros(n, 1), nan(n, 1), ...
      1e-8 + zeros(n, 1), nan(n, 1), VariableNames={ ...
      'tair', 'swd', 'lwd', 'albedo', 'wspd', 'rh', 'psfc', 'ppt', ...
      'rainf', 'snowf'});
   met = icemodel.forcing.helpers.stampMetadata(met);
   met.Properties.UserData = struct( ...
      'source', 'promice', ...
      'precip_policy', ['rainf from PROMICE liquid gauge; snowf and ppt = ' ...
      'NaN placeholders for runtime fill/swap'], ...
      'rainf_source_present', true, ...
      'rainf_observations_present', true, ...
      'rainf_policy', 'rainf observations from liquid gauge');
   met.Properties.UserData.met_resample_policy = "native_15m_unchanged";
   met.Properties.UserData.met_resample_expected_missing_counts = ...
      timetableMissingCounts(met);
   met.Properties.UserData.met_resample_time_semantics = "interval_start";
   met.Properties.UserData.met_resample_support_end_exclusive = ...
      met.Time(end) + minutes(15);
   paths.promice_met = fullfile(met_root, ...
      "met_site_promice_20000101_20000101_15m.mat");
   saveMet(paths.promice_met, met)

   % MERRA carries finite canonical turbulent fluxes and exact MODIS year proof.
   Data = timetable(short_time, -10 + zeros(numel(short_time), 1), ...
      -5 + zeros(numel(short_time), 1), 0.6 + zeros(numel(short_time), 1), ...
      VariableNames={'shf', 'lhf', 'modis'});
   Data = icemodel.forcing.helpers.stampMetadata(Data);
   Data.Properties.UserData = struct( ...
      'merra_flux_sign_convention', 'positive_toward_surface', ...
      'merra_source_time_coordinate', 'native_at_reader', ...
      'merra_time_relabel_policy', ...
      'time_averaged_center_to_interval_start', ...
      'merra_time_upsample_policy', ...
      'zero_order_hold_over_declared_support', ...
      'merra_collection_support_hours', ...
      struct('slv', 1, 'rad', 1, 'flx', 1, 'glc', 3), ...
      'modis_product', 'GEUS Greenland Reflectivity 5km C6', ...
      'modis_status', 'source_coverage', ...
      'modis_coverage_years', 2000);
   paths.merra_data = fullfile(merra_root, ...
      "site_merra2_20000101_20000101.mat");
   saveData(paths.merra_data, Data)

   % MAR includes partial boundary days and complete interior UTC days. The
   % helper produces the same durable hybrid ledger used by real builders.
   mar_time = (datetime(2000, 4, 30, 23, 0, 0, TimeZone="UTC"):hours(1): ...
      datetime(2000, 10, 1, 0, 0, 0, TimeZone="UTC"))';
   groups = findgroups(dateshift(mar_time, 'start', 'day'));
   n_days = max(groups);
   runoff_daily = linspace(0, 0.02, n_days)';
   smb_daily = linspace(-0.01, 0.01, n_days)';
   runoff = runoff_daily(groups) / 24;
   smb = smb_daily(groups) / 24;
   diagnostic_days = linspace(-0.002, 0.002, n_days)';
   subl_evap = diagnostic_days(groups) / 24;
   refreeze_deposition = abs(diagnostic_days(groups)) / 24;
   melt_daily = linspace(0.001, 0.02, n_days)';
   melt = melt_daily(groups) / 24;
   subl = 2e-4 * sin((1:numel(mar_time))' / 24);
   snowd = 0.5 + zeros(numel(mar_time), 1);
   Data = timetable(mar_time, runoff, smb, melt, subl, subl_evap, ...
      refreeze_deposition, snowd, VariableNames={ ...
      'runoff', 'smb', 'melt', 'subl', 'subl_evap', ...
      'refreeze_deposition', 'snowd'});
   Data = icemodel.forcing.helpers.stampMetadata(Data);
   [Data, snow_metadata] = ...
      icemodel.forcing.helpers.applyMarSnowDepthQualityControl(Data);
   [Data, mar_metadata] = ...
      icemodel.forcing.helpers.applyMarDailyQualityControl( ...
      Data, struct('runoff', runoff, 'smb', smb), snow_metadata, sector=1);
   mar_metadata = icemodel.forcing.helpers.marDiagnosticMetadata( ...
      Data, melt, mar_metadata, sector=1);
   Data.Properties.UserData = mar_metadata;
   paths.mar_data = fullfile(mar_root, ...
      "site_mar3.11_20000430_20001001.mat");
   saveData(paths.mar_data, Data)

   % The portable firn manifest references every artifact using production-like
   % source ids and per-leg windows.
   case_entry = struct( ...
      'case_id', 'site', ...
      'case_type', 'firn_observational', ...
      'site_id', 'SITE', ...
      'site_name', 'Synthetic site', ...
      'surface_zone', 'accumulation', ...
      'eval_target', 'surface_temperature', ...
      'permafrost_zone', 'none', ...
      'site_location', struct('lat', 70, 'lon', -40, 'elev', 2000, ...
      'lat_wgs84', 70, 'lon_wgs84', -40), ...
      'period', struct('start', '2000-01-01 00:00:00', ...
      'end', '2000-10-01 00:00:00'), ...
      'evaluation_file', 'site/observations.mat', ...
      'forcing_sources', 'promice', ...
      'eval_sources', {{'promice_obs'; 'mar3.11'; 'merra2'}}, ...
      'comparison_variables', {{'tsfc'; 'runoff'; 'smb'}}, ...
      'observation_variables', struct(), ...
      'colocation', fixtureColocation(), ...
      'native_timestep', 'hourly', ...
      'notes', 'deterministic audit fixture');
   manifest = struct( ...
      'dataset_family', 'promice', ...
      'source_doi', '', ...
      'source_url', 'https://example.invalid', ...
      'source_version', 'fixture', ...
      'retrieval_date', '2000-01-02', ...
      'cases', case_entry);
   paths.manifest = fullfile(family_root, "manifest.json");
   writeJson(paths.manifest, manifest)
end

function [eval_root, input_root, paths] = writeEsmAuditTree(root)
   %WRITEESMAUDITTREE Build one atomic ESM case plus an unrelated met decoy.
   eval_root = fullfile(root, "eval");
   input_root = fullfile(root, "input");
   family_root = fullfile(eval_root, "esm_snowmip");
   case_root = fullfile(family_root, "cdp");
   met_root = fullfile(input_root, "met", "esm_snowmip");
   for folder = [case_root, met_root]
      mkdir(folder)
   end

   % Use one closed 15-minute interval so runtime naming, audit period checks,
   % and source-gap provenance all describe the same exact support.
   times = (datetime(2012, 1, 1, 0, 0, 0, TimeZone="UTC"):minutes(15): ...
      datetime(2012, 1, 2, 23, 45, 0, TimeZone="UTC"))';
   n = numel(times);
   obs = timetable(times, 0.1 + zeros(n, 1), ...
      VariableNames="snow_depth_m");
   obs = icemodel.forcing.helpers.stampMetadata(obs);
   targets = struct('format', 'timeseries', 'data', obs, ...
      'metadata', struct('source_family', 'esm_snowmip', 'station', 'cdp'));
   targets = icemodel.verification.setup.stampArtifactMetadata(targets);
   paths.observations = fullfile(case_root, "observations.mat");
   save(paths.observations, 'targets')

   % The met payload uses every current ESM channel and complete resampling
   % metadata, so the focused test proves full first-class channel inspection.
   met = timetable(times, 263.15 + zeros(n, 1), 100 + zeros(n, 1), ...
      200 + zeros(n, 1), 0.7 + zeros(n, 1), 5 + zeros(n, 1), ...
      80 + zeros(n, 1), 80000 + zeros(n, 1), 1e-8 + zeros(n, 1), ...
      5e-9 + zeros(n, 1), 5e-9 + zeros(n, 1), 0.1 + zeros(n, 1), ...
      VariableNames={'tair', 'swd', 'lwd', 'albedo', 'wspd', 'rh', ...
      'psfc', 'ppt', 'rainf', 'snowf', 'snow_depth'});
   met = icemodel.forcing.helpers.stampMetadata(met);
   met.Properties.UserData = struct( ...
      'source_family', 'esm_snowmip', ...
      'station', 'cdp', ...
      'met_resample_policy', 'native_15m_unchanged', ...
      'met_resample_expected_missing_counts', timetableMissingCounts(met), ...
      'met_resample_time_semantics', 'interval_start', ...
      'met_resample_support_end_exclusive', times(end) + minutes(15));
   met_name = "met_cdp_esm_snowmip_20120101_20120102_15m.mat";
   paths.met_nested = fullfile(met_root, met_name);
   paths.met_flat = fullfile(input_root, "met", met_name);
   paths.decoy = fullfile(met_root, ...
      "met_cdp_esm_snowmip_20100101_20100102_15m.mat");
   saveMet(paths.met_nested, met)
   saveMet(paths.decoy, met)

   % Preserve the intentional atomic schema: observations and forcing share one
   % conversion, so no forcing_sources, colocation, or met_files are declared.
   values = { ...
      'cdp'
      'esm_site'
      'cdp'
      'Col de Porte'
      'land'
      {'seasonal_snow'}
      'none'
      char(fullfile('cdp', 'observations.mat'))
      ''
      '15m'
      struct('start', '2012-01-01 00:00:00', ...
      'end', '2012-01-02 23:45:00')
      {'snow_depth_m'}
      struct('obs_file', char(fullfile('cdp', 'observations.mat')))
      'atomic ESM audit fixture'};
   entry = icemodel.verification.setup.makeCaseManifestEntry(values);
   manifest = icemodel.verification.setup.makeFamilyManifest( ...
      "esm_snowmip", "", "", "fixture", "2012-01-03", entry);
   paths.manifest = fullfile(family_root, "manifest.json");
   icemodel.verification.setup.writeManifest(paths.manifest, manifest);
end

function colocation = fixtureColocation()
   %FIXTURECOLOCATION Return exact file/window declarations for the fixture.
   colocation.promice = struct( ...
      'kind', 'station_met_and_eval', 'staged', true, ...
      'data_files', 'promice/site_promice_20000101_20000101.mat', ...
      'met_files', 'promice/met_site_promice_20000101_20000101_15m.mat', ...
      'window', struct('start', '2000-01-01 00:00:00', ...
      'end', '2000-01-01 03:00:00'));
   colocation.mar = struct( ...
      'kind', 'point_met', 'staged', true, 'source', 'mar', ...
      'source_id', 'mar3.11', 'sample_method', 'nearest', ...
      'data_files', 'mar3.11/site_mar3.11_20000430_20001001.mat', ...
      'window', struct('start', '2000-04-30 23:00:00', ...
      'end', '2000-10-01 00:00:00'));
   colocation.merra = struct( ...
      'kind', 'point_met', 'staged', true, 'source', 'merra', ...
      'source_id', 'merra2', 'sample_method', 'nearest', ...
      'data_files', 'merra2/site_merra2_20000101_20000101.mat', ...
      'window', struct('start', '2000-01-01 00:00:00', ...
      'end', '2000-01-01 03:00:00'));
end

function profile_file = attachMarProfileFixture(input_root, manifest_file)
   %ATTACHMARPROFILEFIXTURE Add one valid MAR RO1 model-output sidecar.
   mar_root = fullfile(input_root, 'userdata', 'mar3.11');
   profile_file = fullfile(mar_root, 'site_mar3.11_density_profiles.mat');
   profile_time = datetime(2000, 5, 1, 12, 0, 0, TimeZone="UTC");
   density = table(repmat("site_20000501", 3, 1), ...
      repmat(profile_time, 3, 1), [0; 1; 2], [300; 500; 700], ...
      repmat("mar3.11", 3, 1), repmat("RO1", 3, 1), ...
      repmat("MARv3.11-ERA5-15km-2000.nc", 3, 1), ...
      repmat("nearest", 3, 1), ...
      'VariableNames', {'profile_id', 'datetime', 'depth', 'density', ...
      'source_id', 'source_variable', 'source_file', 'sample_method'});
   density.Properties.VariableUnits = ...
      {'', '', 'm', 'kg m-3', '', '', '', ''};
   reference = struct('format', 'subsurface_profile_bundle', ...
      'data', struct('density', density), 'metadata', struct());
   save(profile_file, 'reference')

   manifest = jsondecode(fileread(manifest_file));
   manifest.cases.colocation.mar.model_output_files = ...
      'mar3.11/site_mar3.11_density_profiles.mat';
   manifest.cases.colocation.mar.model_output_format = ...
      'subsurface_profile_bundle';
   manifest.cases.colocation.mar.model_output_variables = {'density'};
   writeJson(manifest_file, manifest)
end

function racmo_file = attachRacmoSublFixture(input_root, manifest_file)
   %ATTACHRACMOSUBLFIXTURE Add one legacy manifest-referenced RACMO artifact.
   racmo_root = fullfile(input_root, "userdata", "racmo2.3p3");
   mkdir(racmo_root)
   times = (datetime(2000, 1, 1, 0, 0, 0, TimeZone="UTC"):hours(1): ...
      datetime(2000, 1, 1, 3, 0, 0, TimeZone="UTC"))';
   Data = timetable(times, [-0.002; -0.001; 0; 0.0005], ...
      42 + (1:numel(times))', VariableNames={'subl', 'guard'});
   Data = icemodel.forcing.helpers.stampMetadata(Data, strict=false);
   metadata = struct('sample_method', 'nearest', ...
      'lat_wgs84', 70, 'lon_wgs84', -40);
   Data.Properties.UserData = metadata;
   artifact_metadata = metadata;
   auxiliary = struct('label', "preserve", 'values', [1, 2, 3]);
   racmo_file = string(fullfile(racmo_root, ...
      "site_racmo2.3p3_20000101_20000101.mat"));
   save(racmo_file, 'Data', 'artifact_metadata', 'auxiliary')

   % Add the exact path to the current manifest so source scoping never falls
   % back to filename-token discovery.
   manifest = jsondecode(fileread(manifest_file));
   manifest.cases.colocation.racmo = struct( ...
      'kind', 'point_met', 'staged', true, 'source', 'racmo', ...
      'source_id', 'racmo2.3p3', 'sample_method', 'nearest', ...
      'data_files', 'racmo2.3p3/site_racmo2.3p3_20000101_20000101.mat', ...
      'window', struct('start', '2000-01-01 00:00:00', ...
      'end', '2000-01-01 03:00:00'));
   writeJson(manifest_file, manifest)
end

function met_file = attachMerraMetFixture(input_root, manifest_file, policy)
   %ATTACHMERRAMETFIXTURE Add one manifest-referenced 15-minute MERRA met file.
   met_root = fullfile(input_root, "met", "merra2");
   mkdir(met_root)
   met_time = (datetime(2000, 1, 1, 0, 0, 0, TimeZone="UTC"):minutes(15): ...
      datetime(2000, 1, 1, 3, 45, 0, TimeZone="UTC"))';
   met = timetable(ones(numel(met_time), 1), RowTimes=met_time, ...
      VariableNames="runoff");
   met = icemodel.forcing.helpers.stampMetadata(met, strict=false);
   missing = NaT(0, 1, 'TimeZone', 'UTC');
   metadata = struct( ...
      'merra_source_time_coordinate', 'native_at_reader', ...
      'merra_time_relabel_policy', ...
      'time_averaged_center_to_interval_start', ...
      'merra_time_upsample_policy', ...
      'zero_order_hold_over_declared_support', ...
       'merra_collection_support_hours', ...
       struct('slv', 1, 'rad', 1, 'flx', 1, 'glc', 3), ...
       'merra_tavg3_source_grid_policy', ...
       'native_glc_timestamp_inventory', ...
       'merra_tavg3_expected_source_row_count', 2, ...
       'merra_tavg3_source_row_count', 2, ...
       'merra_tavg3_source_time_gap_count', 0, ...
       'merra_tavg3_missing_source_times', missing, ...
       'met_resample_policy', char(policy), ...
      'met_resample_expected_missing_counts', struct('runoff', 0), ...
      'met_resample_time_semantics', 'interval_start', ...
      'met_resample_support_end_exclusive', met_time(end) + minutes(15), ...
      'sample_method', 'nearest', 'lat_wgs84', 70, 'lon_wgs84', -40);
   met.Properties.UserData = metadata;
   met_file = string(fullfile(met_root, ...
      "met_site_merra2_20000101_20000101_15m.mat"));
   saveMet(met_file, met)

   % Add the exact met path to the existing MERRA leg used for file scoping.
   manifest = jsondecode(fileread(manifest_file));
   manifest.cases.site_location.lat_wgs84 = ...
      manifest.cases.site_location.lat;
   manifest.cases.site_location.lon_wgs84 = ...
      manifest.cases.site_location.lon;
   manifest.cases.colocation.merra.met_files = ...
      'merra2/met_site_merra2_20000101_20000101_15m.mat';
   writeJson(manifest_file, manifest)
end

function met_file = attachMarMetFixture(input_root, manifest_file, data_file)
   %ATTACHMARMETFIXTURE Add derived 15-minute MAR met with copied provenance.
   met_root = fullfile(input_root, "met", "mar3.11");
   mkdir(met_root)
   loaded = load(data_file, 'Data');
   source = loaded.Data(:, ["runoff", "smb"]);
   met_time = (source.Time(1):minutes(15):source.Time(end) + minutes(45))';
   rcm = repelem(source.Variables, 4, 1);
   n = numel(met_time);
   met = timetable(met_time, 260 + zeros(n, 1), 100 + zeros(n, 1), ...
      200 + zeros(n, 1), 0.6 + zeros(n, 1), 5 + zeros(n, 1), ...
      80 + zeros(n, 1), 80000 + zeros(n, 1), zeros(n, 1), ...
      zeros(n, 1), zeros(n, 1), rcm(:, 1), rcm(:, 2), ...
      VariableNames={'tair', 'swd', 'lwd', 'albedo', 'wspd', 'rh', ...
      'psfc', 'ppt', 'rainf', 'snowf', 'runoff', 'smb'});
   met = icemodel.forcing.helpers.stampMetadata(met);

   % Preserve source provenance for traceability while adding the independent
   % derived-met ledger that the audit is expected to enforce.
   metadata = loaded.Data.Properties.UserData;
   metadata.met_resample_policy = 'interval_start_zero_order_hold';
   metadata.met_resample_source_missing_counts = timetableMissingCounts(source);
   metadata.met_resample_expected_missing_counts = timetableMissingCounts(met);
   metadata.met_resample_time_semantics = 'interval_start';
   metadata.met_resample_support_end_exclusive = met.Time(end) + minutes(15);
   metadata.sample_method = 'nearest';
   metadata.lat_wgs84 = 70;
   metadata.lon_wgs84 = -40;
   met.Properties.UserData = metadata;
   met_file = string(fullfile(met_root, ...
      "met_site_mar3.11_20000430_20001001_15m.mat"));
   saveMet(met_file, met)

   % Add the exact met path to the existing MAR leg used for file scoping.
   manifest = jsondecode(fileread(manifest_file));
   manifest.cases.colocation.mar.met_files = ...
      'mar3.11/met_site_mar3.11_20000430_20001001_15m.mat';
   writeJson(manifest_file, manifest)
end

function saveData(pathname, Data)
   %SAVEDATA Save one metadata-synchronized userdata artifact.
   artifact_metadata = Data.Properties.UserData;
   save(pathname, 'Data', 'artifact_metadata')
end

function saveMet(pathname, met)
   %SAVEMET Save one metadata-synchronized met artifact.
   artifact_metadata = met.Properties.UserData;
   save(pathname, 'met', 'artifact_metadata')
end

function counts = timetableMissingCounts(T)
   %TIMETABLEMISSINGCOUNTS Mirror writer provenance for the compact fixture.
   counts = struct();
   for name = reshape(string(T.Properties.VariableNames), 1, [])
      values = T.(char(name));
      if isnumeric(values)
         counts.(char(name)) = nnz(~isfinite(values));
      end
   end
end

function writeTinyMerraGlcSource(root, day)
   %WRITETINYMERRAGLCSOURCE Create one official-shape daily native tavg3 file.
   folder = fullfile(root, 'glc');
   mkdir(folder)
   filename = fullfile(folder, ...
      "MERRA2_400.tavg3_2d_glc_Nx." + string(day, 'yyyyMMdd') + ".nc4");
   nccreate(filename, 'RUNOFF', ...
      'Dimensions', {'lon', 1, 'lat', 1, 'time', 8});
   ncwrite(filename, 'RUNOFF', reshape(1:8, 1, 1, []));
   ncwriteatt(filename, 'RUNOFF', 'units', 'kg m-2 s-1');
   nccreate(filename, 'time', 'Dimensions', {'time', 8});
   ncwrite(filename, 'time', (0:180:1260)');
   ncwriteatt(filename, 'time', 'units', ...
      "minutes since " + string(day + minutes(90), 'yyyy-MM-dd HH:mm:ss'));
end

function writeLegacyMerraMet(filename, proven, include_merra)
   %WRITELEGACYMERRAMET Create an hourly-ramped legacy 15-minute MERRA artifact.
   if nargin < 3
      include_merra = true;
   end
   source_time = (datetime(2000, 1, 1, TimeZone="UTC"):hours(1): ...
      datetime(2000, 1, 1, 3, 0, 0, TimeZone="UTC"))';
   source = timetable((1:4)', RowTimes=source_time, ...
      VariableNames="runoff");
   legacy_time = (source_time(1):minutes(15):source_time(end))';
   met = retime(source, legacy_time, 'linear');
   metadata = struct( ...
      'met_resample_policy', 'linear_adjacent_finite_only', ...
      'met_resample_source_row_count', height(source), ...
      'met_resample_source_cadence_seconds', 3600, ...
      'met_resample_source_time_gap_count', 0, ...
      'met_resample_expected_missing_counts', struct('runoff', 0));
   if include_merra
      % Current MERRA timing markers accompany the optional native-grid proof.
      metadata.merra_source_time_coordinate = 'native_at_reader';
      metadata.merra_time_relabel_policy = ...
         'time_averaged_center_to_interval_start';
      metadata.merra_time_upsample_policy = ...
         'zero_order_hold_over_declared_support';
      metadata.merra_collection_support_hours = ...
         struct('slv', 1, 'rad', 1, 'flx', 1, 'glc', 3);
   end
   if proven
      % The two saved 00/03 UTC rows are independently present in native glc.
      metadata.merra_tavg3_source_grid_policy = ...
         'native_glc_timestamp_inventory';
      metadata.merra_tavg3_expected_source_row_count = 2;
      metadata.merra_tavg3_source_row_count = 2;
      metadata.merra_tavg3_source_time_gap_count = 0;
      metadata.merra_tavg3_missing_source_times = ...
         NaT(0, 1, 'TimeZone', 'UTC');
   end
   met.Properties.UserData = metadata;
   save(filename, 'met')
end

function writeJson(pathname, value)
   %WRITEJSON Replace one compact fixture manifest deterministically.
   fid = fopen(pathname, 'w');
   if fid < 0
      error('test:verificationArtifactAudit:manifestOpen', ...
         'could not open fixture manifest: %s', pathname)
   end
   cleanup = onCleanup(@() fclose(fid));
   fwrite(fid, jsonencode(value, PrettyPrint=true), 'char');
   clear cleanup
end

function snapshot = fileSnapshot(roots)
   %FILESNAPSHOT Record path, bytes, and mtime for non-mutation proof.
   rows = cell(numel(roots), 1);
   for k = 1:numel(roots)
      entries = dir(fullfile(roots(k), '**', '*'));
      entries = entries(~[entries.isdir]);
      rows{k} = table(string(fullfile({entries.folder}, {entries.name}))', ...
         [entries.bytes]', [entries.datenum]', VariableNames={ ...
         'path', 'bytes', 'datenum'});
   end
   snapshot = sortrows(vertcat(rows{:}), 'path');
end

function bytes = fileBytes(filename)
   %FILEBYTES Return exact bytes so a dry-run mutation cannot hide in metadata.
   fid = fopen(filename, 'r');
   assert(fid >= 0, 'could not open fixture artifact')
   cleanup = onCleanup(@() fclose(fid));
   bytes = fread(fid, Inf, '*uint8');
end

function record = channelRecord(channels, source, channel)
   %CHANNELRECORD Return one source/channel diagnostic from the structured report.
   keep = string({channels.source}) == source ...
      & string({channels.channel}) == channel;
   record = channels(find(keep, 1));
end

function verifyCodes(testCase, actual, expected)
   %VERIFYCODES Assert each independent defect class appears in one audit pass.
   for code = expected
      testCase.verifyTrue(any(actual == code), code);
   end
end
