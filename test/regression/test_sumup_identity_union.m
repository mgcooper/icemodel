function tests = test_sumup_identity_union
   %TEST_SUMUP_IDENTITY_UNION Prove the bounded canonical SUMup identity union.
   %
   % This source-backed restage is intentionally opt-in because it rereads the
   % multi-million-row SUMup files for every retained default case. Run it once
   % with ICEMODEL_RUN_SUMUP_IDENTITY_UNION=1 before canonical replacement.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % Install project dependencies for the explicitly enabled bounded restage.
   [~, ~, ~, ~, cleanup] = ...
      icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;
   testCase.TestData.repo_root = string(icemodel.internal.fullpath());
   testCase.TestData.cache = fullfile(testCase.TestData.repo_root, ...
      'data', 'verification', 'sumup');
end

function teardownOnce(testCase)
   % Release the shared test environment after the opt-in audit procedure.
   testCase.TestData.cleanup = [];
end

function test_default_restage_preserves_union_and_every_case_contributes(testCase)
   % The canonical 47-case SUMup set plus external Humphrey must preserve every
   % audited identity, report exact per-artifact removals, and pass artifact QA.
   testCase.assumeEqual(string(getenv( ...
      'ICEMODEL_RUN_SUMUP_IDENTITY_UNION')), "1", ...
      ['Set ICEMODEL_RUN_SUMUP_IDENTITY_UNION=1 for the bounded, ' ...
      'fixture-heavy SUMup restage audit.']);
   verifySourceHashes(testCase, testCase.TestData.cache);

   tmp = string(tempname);
   mkdir(tmp)
   testCase.addTeardown(@() rmdir(tmp, 's'));
   eval_root = fullfile(tmp, 'eval');
   input_root = fullfile(tmp, 'input');
   copyAnchorManifests(testCase.TestData.repo_root, eval_root);

   manifest = icemodel.verification.setup.importSumup( ...
      testCase.TestData.cache, evaluation_data_root=eval_root, ...
      input_data_root=input_root, build_forcing=false, ...
      overwrite_family=true, skip_missing=false);
   case_ids = string({manifest.cases.case_id});
   testCase.verifyEqual(numel(case_ids), 47);
   testCase.verifyFalse(any(ismember(["serb", "tasu", "thul", "thul2", ...
      "zaca", "zacl", "humphrey", "s23", "lynl", "lynt", "nukb"], ...
      case_ids)));
   testCase.verifyTrue(all(ismember(["dy2", "dye-2long"], case_ids)));

   bundles = ["density", "subsurface_temperature", "smb"];
   source_variables = ["density", "temperature", "SMB"];
   expected_unique = [913771, 1717279, 55939];
   expected_artifact_removed = [171, 4628, 360];
   parts = struct('density', {{}}, 'subsurface_temperature', {{}}, ...
      'smb', {{}});
   owners = parts;
   artifact_removed = zeros(1, numel(bundles));

   % Load only the newly staged observation artifacts and retain numeric identity
   % projections for one cross-case union pass.
   for k = 1:numel(manifest.cases)
      loaded = load(fullfile(eval_root, 'sumup', ...
         manifest.cases(k).evaluation_file), 'targets');
      for b = 1:numel(bundles)
         bundle = bundles(b);
         record = loaded.targets.data.(bundle);
         meta = loaded.targets.metadata;
         artifact_removed(b) = artifact_removed(b) + ...
            meta.(char(bundle + "_duplicate_rows_removed"));
         if isempty(record)
            continue
         end
         identity = identityProjection(record, bundle);
         parts.(bundle){end + 1} = identity;
         owners.(bundle){end + 1} = repmat(k, height(identity), 1);
      end
   end
   testCase.verifyEqual(artifact_removed, expected_artifact_removed);

   % Humphrey is externally owned but must participate in the preservation and
   % exclusivity comparison so the omitted SUMup alias cannot hide lost rows.
   site = icemodel.verification.setup.researchSiteCatalog("humphrey");
   [external, external_meta] = ...
      icemodel.verification.setup.buildSumupObservations( ...
      [site.lat_wgs84, site.lon_wgs84], ...
      source_dir=testCase.TestData.cache);
   testCase.verifyEqual( ...
      external_meta.subsurface_temperature_duplicate_rows_removed, 119);
   for b = 1:numel(bundles)
      bundle = bundles(b);
      if isempty(external.(bundle))
         continue
      end
      identity = identityProjection(external.(bundle), bundle);
      parts.(bundle){end + 1} = identity;
      owners.(bundle){end + 1} = zeros(height(identity), 1);
   end

   % Group normalized identities once per variable. A group owned by exactly one
   % positive case index is that retained SUMup case's exclusive contribution.
   contributes = false(numel(manifest.cases), 1);
   for b = 1:numel(bundles)
      bundle = bundles(b);
      identity = vertcat(parts.(bundle){:});
      owner = vertcat(owners.(bundle){:});
      comparison = normalizeMissingIdentity(identity);
      [~, ~, group] = unique(comparison, 'rows');
      testCase.verifyEqual(max(group), expected_unique(b));
      group_min = accumarray(group, owner, [], @min);
      group_max = accumarray(group, owner, [], @max);
      exclusive = unique(group_min(group_min == group_max & group_min > 0));
      contributes(exclusive) = true;

      % The production helper independently agrees that the union contains the
      % exact expected scientific identities for this source variable.
      [~, ~, unique_rows] = ...
         icemodel.verification.setup.deduplicateSumupRecords( ...
         identity, source_variables(b));
      testCase.verifyEqual(unique_rows, expected_unique(b));
   end
   testCase.verifyTrue(all(contributes), ...
      "Every retained canonical SUMup case must add at least one identity.");

   report = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families="sumup", report_dir="");
   testCase.verifyTrue(report.summary.passed, ...
      strjoin(string({report.findings.code}), newline));
end

function verifySourceHashes(testCase, cache)
   %VERIFYSOURCEHASHES Pin the immutable release bytes used by identity counts.
   files = ["SUMup_2025_density_greenland.nc", ...
      "SUMup_2025_temperature_greenland.nc", ...
      "SUMup_2025_SMB_greenland.nc"];
   expected = [ ...
      "294674c459bdcec99e0010901e644d634777431e61c741554a2df8de81af3ba4"
      "1eb0510de5a72cefa943840eaf282c897d82db06e67a3e3582ae050c97c2c6d6"
      "72de0bbb8fb7ea0519cf56619d5dba17e6544f37013d5c80241a33629f757562"];
   for k = 1:numel(files)
      pathname = fullfile(cache, files(k));
      testCase.assumeTrue(isfile(pathname), ...
         sprintf('SUMup source absent: %s', pathname));
      testCase.verifyEqual( ...
         icemodel.verification.setup.fileSha256(pathname), expected(k));
   end
end

function copyAnchorManifests(repo_root, eval_root)
   %COPYANCHORMANIFESTS Seed only the read-only mixed catalog into a temp root.
   families = setdiff( ...
      icemodel.verification.namelists.firndatasetfamily(), "sumup", 'stable');
   for family = reshape(families, 1, [])
      source = fullfile(repo_root, 'data', 'eval', family, 'manifest.json');
      if ~isfile(source)
         continue
      end
      destination = fullfile(eval_root, family);
      mkdir(destination)
      copyfile(source, fullfile(destination, 'manifest.json'));
   end
end

function identity = identityProjection(record, bundle)
   %IDENTITYPROJECTION Reverse datetime shaping and strip non-identity metadata.
   epoch = datetime(1900, 1, 1, 'TimeZone', 'UTC');
   switch bundle
      case "density"
         record.timestamp = days(record.datetime - epoch);
         record.datetime = [];
         names = ["latitude", "longitude", "start_depth", "stop_depth", ...
            "midpoint", "timestamp", "reference_key", "method_key", ...
            "density", "error"];
      case "subsurface_temperature"
         record = timetable2table(record, 'ConvertRowTimes', true);
         time_name = string(record.Properties.VariableNames{1});
         record.timestamp = days(record.(time_name) - epoch);
         record.(time_name) = [];
         record = renamevars(record, "subsurface_temperature", "temperature");
         names = ["latitude", "longitude", "depth", "duration", ...
            "timestamp", "reference_key", "method_key", "temperature", ...
            "error"];
      case "smb"
         record.start_date = days(record.start_date - epoch);
         record.end_date = days(record.end_date - epoch);
         names = ["latitude", "longitude", "start_date", "end_date", ...
            "start_year", "end_year", "reference_key", "method_key", ...
            "smb", "error"];
   end

   % Rebuild a plain numeric table so vertical concatenation is independent of
   % per-artifact custom table metadata.
   selected = record(:, cellstr(names));
   identity = array2table(selected{:,:}, ...
      'VariableNames', selected.Properties.VariableNames);
end

function identity = normalizeMissingIdentity(identity)
   %NORMALIZEMISSINGIDENTITY Make pre-R2026a UNIQUE grouping missing-aware.
   names = string(identity.Properties.VariableNames);
   for name = names
      values = identity.(name);
      missing = ismissing(values);
      if ~any(missing)
         continue
      end
      identity.(name + "_is_missing") = missing;
      values(missing) = 0;
      identity.(name) = values;
   end
end
