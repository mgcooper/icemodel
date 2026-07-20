function repairOmittedRacmoPrecipitationArtifacts()
   %REPAIROMITTEDRACMOPRECIPITATIONARTIFACTS Repair five missed artifacts.
   %
   % Applies the existing RACMO precipitation/QC contract to the exact RetMIP
   % and IMAU manifest references omitted from the original PROMICE-scoped
   % write pass. Sublimation values and markers are deliberately untouched.

   repo_root = string(icemodel.internal.fullpath());
   data_root = fullfile(repo_root, "data");
   input_root = fullfile(data_root, "input");
   eval_root = fullfile(data_root, "eval");
   qa_root = fullfile(data_root, "preview", "qa");
   evidence_root = fullfile(repo_root, "sandbox", "verification", ...
      "icemodel-hfc.2.38", "evidence");
   icemodel.helpers.ensureDirExists(qa_root)
   icemodel.helpers.ensureDirExists(evidence_root)

   families = ["retmip", "imau"];
   manifests = fullfile(eval_root, families, "manifest.json");
   manifest_hashes = strings(size(manifests));
   for k = 1:numel(manifests)
      manifest_hashes(k) = ...
         icemodel.verification.setup.fileSha256(manifests(k));
   end

   % Persist the complete dry plan before mutation so a failed write remains
   % diagnosable. Disabling the later sign repair keeps this transaction narrow.
   dry = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family=families, ...
      source_id="racmo2.3p3", repair_racmo_subl=false);
   assertExactFiles(dry)
   assertPrecipitationOnly(dry)
   preflight = repairRows(dry, "dry");
   writetable(preflight, fullfile(evidence_root, ...
      "racmo_precipitation_preflight.csv"))
   before = payloadSnapshots(dry);

   write = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family=families, ...
      source_id="racmo2.3p3", repair_racmo_subl=false, dry_run=false);
   second = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family=families, ...
      source_id="racmo2.3p3", repair_racmo_subl=false, dry_run=false);
   assertExactFiles(write)
   assert(all(string({write.records.status}) == "repaired"), ...
      'the five omitted RACMO artifacts were not all repaired')
   assert(all(string({second.records.status}) == "unchanged"), ...
      'RACMO precipitation repair pass two was not unchanged')
   assertActionEquality(dry, write)
   validatePayloads(write, second, before)

   % Manifests and exact references are inputs to the repair, never outputs.
   for k = 1:numel(manifests)
      assert(icemodel.verification.setup.fileSha256(manifests(k)) ...
         == manifest_hashes(k), 'RACMO repair changed a family manifest')
   end
   audit = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families=families);
   source = string({audit.findings.source});
   channel = string({audit.findings.channel});
   ppt_finding = ismember(lower(source), ["racmo", "racmo2.3p3"]) ...
      & ismember(channel, ["ppt", "precip"]);
   assert(~any(ppt_finding), ...
      'focused audit still reports a RACMO precipitation finding')

   evidence = [repairRows(dry, "dry"); repairRows(write, "write"); ...
      repairRows(second, "second")];
   writetable(evidence, fullfile(qa_root, ...
      "racmo_omitted_precipitation_repair.csv"))
   fprintf('Repaired five omitted RetMIP/IMAU RACMO precipitation artifacts.\n');
end

function assertExactFiles(report)
   %ASSERTEXACTFILES Require only the five reviewed manifest references.
   expected = sort([ ...
      "dye2_long_racmo2.3p3_20120101_20150502.mat"
      "dye2_2016_racmo2.3p3_20160502_20161028.mat"
      "summit_racmo2.3p3_20120101_20150308.mat"
      "fa_racmo2.3p3_20140412_20141202.mat"
      "s22_racmo2.3p3_20160814_20181231.mat"]);
   actual = strings(numel(report.records), 1);
   for k = 1:numel(report.records)
      [~, name, ext] = fileparts(report.records(k).filename);
      actual(k) = string(name) + string(ext);
   end
   assert(isequal(sort(actual), expected), ...
      'manifest-derived RACMO precipitation inventory is not the exact five')
end

function assertPrecipitationOnly(report)
   %ASSERTPRECIPITATIONONLY Reject sign or unrelated dry-run transformations.
   allowed_actions = ["rename_racmo_precip_to_ppt", ...
      "zero_negative_racmo_ppt", "stamp_racmo_ppt_qc", ...
      "restamp_metadata", "sync_artifact_metadata"];
   for k = 1:numel(report.records)
      record = report.records(k);
      assert(string(record.status) == "would_repair", ...
         'omitted RACMO precipitation artifact is not repairable')
      assert(all(ismember(string(record.actions), allowed_actions)), ...
         'dry plan contains a non-precipitation action')
      assert(all(ismember(string(record.changed_variables), ...
         ["precip", "ppt"])), 'dry plan changes a non-precipitation variable')
      metadata = string(record.changed_metadata_fields);
      assert(all(startsWith(metadata, "racmo_ppt_qc_")), ...
         'dry plan changes unrelated metadata')
      assert(~any(contains(string(record.actions), "racmo_subl")), ...
         'precipitation pre-restamp planned a sublimation action')
   end
end

function snapshots = payloadSnapshots(report)
   %PAYLOADSNAPSHOTS Retain exact pre-write Time, subl, and unrelated values.
   snapshots = cell(numel(report.records), 1);
   for k = 1:numel(report.records)
      loaded = load(report.records(k).filename, 'Data', 'artifact_metadata');
      snapshots{k} = loaded;
   end
end

function assertActionEquality(dry, write)
   %ASSERTACTIONEQUALITY Require the dry and write plans to be identical.
   for k = 1:numel(dry.records)
      assert(isequal(string(dry.records(k).actions), ...
         string(write.records(k).actions)), ...
         'RACMO precipitation dry/write actions differ')
      assert(isequal(string(dry.records(k).changed_variables), ...
         string(write.records(k).changed_variables)), ...
         'RACMO precipitation dry/write variable ledgers differ')
      assert(isequal(string(dry.records(k).changed_metadata_fields), ...
         string(write.records(k).changed_metadata_fields)), ...
         'RACMO precipitation dry/write metadata ledgers differ')
   end
end

function validatePayloads(write, second, before)
   %VALIDATEPAYLOADS Prove precipitation-only change and pass-two hashes.
   required = ["racmo_ppt_qc_method", "racmo_ppt_qc_stage", ...
      "racmo_ppt_qc_source_variable", "racmo_ppt_qc_basis", ...
      "racmo_ppt_qc_status", "racmo_ppt_qc_replaced_count", ...
      "racmo_ppt_qc_input_minimum", "racmo_ppt_qc_output_minimum"];
   for k = 1:numel(write.records)
      after = load(write.records(k).filename, 'Data', 'artifact_metadata');
      prior = before{k};
      assert(isequaln(after.Data.Properties.RowTimes, ...
         prior.Data.Properties.RowTimes), 'RACMO repair changed Time')
      assert(isequaln(after.Data.subl, prior.Data.subl), ...
         'RACMO precipitation repair changed subl')
      assert(~isfield(after.artifact_metadata, ...
         'racmo_subl_sign_convention') ...
         && ~isfield(after.artifact_metadata, ...
         'racmo_subl_native_sign_convention'), ...
         'RACMO precipitation repair added sign markers')
      assert(ismember("ppt", string(after.Data.Properties.VariableNames)) ...
         && ~ismember("precip", string(after.Data.Properties.VariableNames)), ...
         'RACMO precipitation was not canonicalized to ppt')
      assert(all(after.Data.ppt(isfinite(after.Data.ppt)) >= 0), ...
         'canonical RACMO ppt remains negative')
      assert(all(isfield(after.artifact_metadata, required)), ...
         'RACMO precipitation QC metadata is incomplete')

      % The repair may add only the reviewed precipitation provenance fields;
      % every pre-existing unrelated top-level metadata value remains exact.
      metadata_fields = union(string(fieldnames(prior.artifact_metadata)), ...
         string(fieldnames(after.artifact_metadata)), "stable");
      unrelated_metadata = metadata_fields( ...
         ~startsWith(metadata_fields, "racmo_ppt_qc_"));
      for name = reshape(unrelated_metadata, 1, [])
         assert(isfield(prior.artifact_metadata, name) ...
            && isfield(after.artifact_metadata, name) ...
            && isequaln(prior.artifact_metadata.(name), ...
            after.artifact_metadata.(name)), ...
            'RACMO repair changed unrelated metadata %s', name)
      end

      % Every common non-precipitation variable is exact; the helper separately
      % rejects any unexpected added/removed variable name.
      prior_names = string(prior.Data.Properties.VariableNames);
      after_names = string(after.Data.Properties.VariableNames);
      common = setdiff(intersect(prior_names, after_names, "stable"), ...
         ["precip", "ppt"]);
      for name = reshape(common, 1, [])
         assert(isequaln(prior.Data.(name), after.Data.(name)), ...
            'RACMO repair changed unrelated variable %s', name)
      end
      assert(write.records(k).unrelated_payload_preserved, ...
         'shared repair did not preserve unrelated payload')
      assert(second.records(k).hash_before == write.records(k).hash_after ...
         && second.records(k).hash_after == write.records(k).hash_after, ...
         'RACMO precipitation pass-two hash changed')
   end
end

function rows = repairRows(report, phase)
   %REPAIRROWS Flatten one exact repair phase into durable evidence.
   n = numel(report.records);
   phase = repmat(string(phase), n, 1);
   artifact_path = strings(n, 1);
   status = strings(n, 1);
   actions = strings(n, 1);
   changed_variables = strings(n, 1);
   changed_metadata_fields = strings(n, 1);
   hash_before = strings(n, 1);
   hash_after = strings(n, 1);
   for k = 1:n
      [~, name, ext] = fileparts(report.records(k).filename);
      artifact_path(k) = "data/input/userdata/racmo2.3p3/" ...
         + string(name) + string(ext);
      status(k) = string(report.records(k).status);
      actions(k) = strjoin(string(report.records(k).actions), ";");
      changed_variables(k) = strjoin(string( ...
         report.records(k).changed_variables), ";");
      changed_metadata_fields(k) = strjoin(string( ...
         report.records(k).changed_metadata_fields), ";");
      hash_before(k) = string(report.records(k).hash_before);
      hash_after(k) = string(report.records(k).hash_after);
   end
   rows = table(phase, artifact_path, status, actions, changed_variables, ...
      changed_metadata_fields, hash_before, hash_after);
end
