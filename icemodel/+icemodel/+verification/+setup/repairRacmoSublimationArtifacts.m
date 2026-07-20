function repairRacmoSublimationArtifacts(candidate_root)
   %REPAIRRACMOSUBLIMATIONARTIFACTS Repair referenced RACMO artifacts.
   %
   % Repairs only exact manifest-referenced RACMO userdata in the isolated
   % seasonal candidate and canonical repository roots. Compact pre/post hash,
   % value-change, idempotence, affected-figure, and focused-audit evidence is
   % written under data/preview/qa.

   arguments
      candidate_root (1, 1) string
   end

   repo_root = string(icemodel.internal.fullpath());
   data_root = fullfile(repo_root, "data");
   code_root = fullfile(repo_root, "icemodel");
   qa_root = fullfile(data_root, "preview", "qa");
   evidence_root = fullfile(repo_root, "sandbox", "verification", ...
      "icemodel-hfc.2.38", "evidence");
   addpath(code_root)
   icemodel.helpers.ensureDirExists(qa_root)
   icemodel.helpers.ensureDirExists(evidence_root)

   % A nonempty family list activates exact manifest-file scoping in the shared
   % repair helper. Families without RACMO references contribute no records.
   families = ["promice", "retmip", "imau", "sumup", "research_site"];
   roots = struct( ...
      'kind', {"candidate", "canonical"}, ...
      'input', {fullfile(candidate_root, "input"), ...
      fullfile(data_root, "input")}, ...
      'eval', {fullfile(candidate_root, "eval"), ...
      fullfile(data_root, "eval")});
   for k = 1:numel(roots)
      assert(isfolder(roots(k).input) && isfolder(roots(k).eval), ...
         'repair root is incomplete: %s', roots(k).kind)
   end

   family_sets = {"promice"; families};

   % Inventory both physical roots before the first mutation. A legacy artifact
   % has neither marker and plans the one-time flip; paired current markers are
   % unchanged; partial or unknown markers fail this preflight.
   dry_reports = cell(numel(roots), 1);
   preflight = cell(numel(roots), 1);
   for k = 1:numel(roots)
      dry_reports{k} = ...
         icemodel.verification.setup.repairRcmArtifactMetadata( ...
         roots(k).input, eval_root=roots(k).eval, ...
         dataset_family=family_sets{k}, source_id="racmo2.3p3");
      assert(all(ismember(string({dry_reports{k}.records.status}), ...
         ["would_repair", "unchanged"])), ...
         'RACMO dry run found an ambiguous or unrepairable artifact')
      assertOnlySublimationChanges(dry_reports{k})
      preflight{k} = preflightRows(roots(k), dry_reports{k});
   end
   preflight = vertcat(preflight{:});
   assert(height(preflight(preflight.root_kind == "candidate", :)) == 42 ...
      && height(preflight(preflight.root_kind == "canonical", :)) == 47, ...
      'preflight does not contain all 89 exact manifest references')
   writetable(preflight, fullfile(evidence_root, ...
      "racmo_repair_preflight.csv"))
   fprintf('RACMO preflight: %d canonical, %d legacy-native, 0 ambiguous.\n', ...
      nnz(preflight.sign_state == "canonical"), ...
      nnz(preflight.sign_state == "legacy_native_unmarked"));

   % Pass one repairs only remaining legacy artifacts; pass two is the
   % byte-stable idempotence proof. The shared helper verifies unrelated payload.
   inventory = cell(numel(roots), 1);
   audit_rows = cell(numel(roots), 1);
   for k = 1:numel(roots)
      first = icemodel.verification.setup.repairRcmArtifactMetadata( ...
         roots(k).input, eval_root=roots(k).eval, ...
         dataset_family=family_sets{k}, source_id="racmo2.3p3", ...
         dry_run=false);
      second = icemodel.verification.setup.repairRcmArtifactMetadata( ...
         roots(k).input, eval_root=roots(k).eval, ...
         dataset_family=family_sets{k}, source_id="racmo2.3p3", ...
         dry_run=false);
      assert(~isempty(first.records), 'no manifest-referenced RACMO artifacts')
      assert(all(ismember(string({first.records.status}), ...
         ["repaired", "unchanged"])), 'RACMO repair did not converge')
      assert(all(string({second.records.status}) == "unchanged"), ...
         'RACMO repair pass two was not unchanged')
      assert(all([first.records.unrelated_payload_preserved]), ...
         'RACMO repair changed unrelated payload')
      inventory{k} = repairRows(roots(k), first, second);

      % Source-aware QA must no longer report either missing or noncanonical
      % RACMO sublimation markers; unrelated findings remain outside this Bead.
      audit = icemodel.verification.auditArtifacts( ...
         evaluation_data_root=roots(k).eval, ...
         input_data_root=roots(k).input, families=family_sets{k});
      codes = string({audit.findings.code});
      sign_findings = startsWith(codes, "racmo_subl_sign");
      assert(~any(sign_findings), ...
         'focused artifact QA still reports a RACMO sublimation sign finding')
      audit_rows{k} = table(roots(k).kind, numel(codes), nnz(sign_findings), ...
         VariableNames=["root_kind", "finding_count", ...
         "racmo_subl_sign_finding_count"]);
   end

   inventory = vertcat(inventory{:});
   audit_summary = vertcat(audit_rows{:});
   figure_inventory = affectedFigureRows(candidate_root);
   flux_sign_evidence = fluxSignEvidence(candidate_root);
   assert(height(inventory(inventory.root_kind == "candidate", :)) == 42, ...
      'candidate manifest no longer references the reviewed 42 RACMO artifacts')
   assert(height(inventory(inventory.root_kind == "canonical", :)) == 47, ...
      'canonical manifests no longer reference the reviewed 47 RACMO artifacts')
   assert(height(figure_inventory) == 42, ...
      'candidate figure inventory no longer contains 42 RACMO subl groups')
   inventory_file = fullfile(qa_root, ...
      "racmo_sublimation_repair_inventory.csv");
   if isfile(inventory_file) && all(inventory.first_status == "unchanged")
      % Preserve first-pass pre/post hashes on later idempotence replays while
      % proving every current hash still equals the recorded repaired hash.
      prior = readtable(inventory_file, TextType="string");
      [prior, current] = alignEvidenceRows(prior, inventory);
      assert(all(prior.hash_after == current.hash_after), ...
         'current RACMO hashes differ from preserved first-pass evidence')
      inventory = prior;
   else
      writetable(inventory, inventory_file)
   end
   writetable(figure_inventory, fullfile(qa_root, ...
      "racmo_sublimation_figure_inventory.csv"))
   writetable(audit_summary, fullfile(qa_root, ...
      "racmo_sublimation_audit_summary.csv"))
   writetable(flux_sign_evidence, fullfile(qa_root, ...
      "verification_flux_sign_evidence.csv"))
   fprintf('RACMO sublimation repair complete: %d artifacts, %d figures.\n', ...
      height(inventory), height(figure_inventory));
end

function rows = preflightRows(root, report)
   %PREFLIGHTROWS Persist exact pre-mutation hashes and inferred sign state.
   n = numel(report.records);
   root_kind = repmat(string(root.kind), n, 1);
   artifact_path = strings(n, 1);
   hash_before = strings(n, 1);
   sign_state = strings(n, 1);
   planned_actions = strings(n, 1);
   for k = 1:n
      artifact_path(k) = root_kind(k) + "/" + string( ...
         erase(report.records(k).filename, string(root.input) + filesep));
      hash_before(k) = string(report.records(k).hash_before);
      planned_actions(k) = strjoin(string(report.records(k).actions), ";");
      if string(report.records(k).status) == "unchanged"
         sign_state(k) = "canonical";
      else
         sign_state(k) = "legacy_native_unmarked";
      end
   end
   rows = table(root_kind, artifact_path, hash_before, sign_state, ...
      planned_actions);
end

function [prior, current] = alignEvidenceRows(prior, current)
   %ALIGNEVIDENCEROWS Sort two repair ledgers by their exact artifact key.
   [~, prior_order] = sort(prior.root_kind + "/" + prior.artifact_path);
   [~, current_order] = sort(current.root_kind + "/" + current.artifact_path);
   prior = prior(prior_order, :);
   current = current(current_order, :);
   assert(isequal(prior.root_kind, current.root_kind) ...
      && isequal(prior.artifact_path, current.artifact_path), ...
      'current RACMO inventory differs from preserved first-pass evidence')
end

function evidence = fluxSignEvidence(candidate_root)
   %FLUXSIGNEVIDENCE Formalize all-source heat and sublimation sign decisions.
   filename = fullfile(candidate_root, "preview", "qa", ...
      "flux_sign_pair_stats.csv");
   assert(isfile(filename), 'candidate flux-sign pair evidence is missing')
   source = readtable(filename, TextType="string");
   % Require the complete reviewed matrix so missing source/site evidence cannot
   % silently inherit a favorable sign decision from the remaining rows.
   sites = ["cen", "kanl", "kanm"];
   heat_sources = ["mar3_11", "merra2", "racmo2_3p3"];
   for site = sites
      for model = heat_sources
         for channel = ["shf", "lhf"]
            assert(nnz(source.site == site & source.source == model ...
               & source.channel == channel) == 1, ...
               'all-source heat-flux evidence matrix is incomplete')
         end
      end
      assert(nnz(source.site == site ...
         & source.source == "racmo2_3p3_vs_mar" ...
         & source.channel == "subl") == 1, ...
         'MAR/RACMO sublimation evidence matrix is incomplete')
   end
   assert(height(source) == 21, ...
      'flux-sign evidence contains an unexpected row')
   evidence = source;
   names = string(evidence.Properties.VariableNames);
   names(names == "correlation_current") = "correlation_pre_repair";
   names(names == "correlation_flipped") = ...
      "alternative_negated_correlation";
   evidence.Properties.VariableNames = cellstr(names);
   evidence.reference_source = repmat("promice", height(evidence), 1);
   evidence.correlation_post_repair = evidence.correlation_pre_repair;
   evidence.pre_orientation = repmat( ...
      "positive_toward_surface", height(evidence), 1);
   evidence.post_orientation = evidence.pre_orientation;
   evidence.decision = repmat("canonical_supported_no_mutation", ...
      height(evidence), 1);
   evidence.evidence_basis = repmat( ...
      "daily_mean_against_promice_observations", height(evidence), 1);

   % Weak CEN correlations do not identify a sign error and must never trigger
   % mutation; KAN_L/KAN_M provide the discriminating positive-orientation rows.
   heat = ismember(evidence.channel, ["shf", "lhf"]);
   ambiguous = heat & evidence.site == "cen";
   evidence.decision(ambiguous) = ...
      "ambiguous_weak_model_difference_no_mutation";

   % The old table records the tested hypothetical RACMO negation. Recompute
   % against the now-repaired values so the durable post column is independent
   % evidence rather than a relabeled prediction.
   subl = evidence.channel == "subl";
   evidence.reference_source(subl) = "mar3.11";
   evidence.pre_orientation(subl) = ...
      "negative_loss_positive_deposition";
   evidence.post_orientation(subl) = ...
      "positive_loss_negative_deposition";
   evidence.evidence_basis(subl) = "daily_mean_racmo_against_mar";
   evidence.decision(subl) = "repaired_to_canonical_positive_loss";
   for row = find(subl)'
      post = postRepairSublimationCorrelation(candidate_root, ...
         evidence.site(row));
      assert(evidence.correlation_pre_repair(row) < 0 && post > 0, ...
         'MAR/RACMO sublimation did not reverse to positive correlation')
      assert(abs(post - evidence.alternative_negated_correlation(row)) < 1e-10, ...
         'post-repair sublimation differs from the source-backed prediction')
      evidence.correlation_post_repair(row) = post;
      evidence.alternative_negated_correlation(row) = -post;
   end
end

function value = postRepairSublimationCorrelation(candidate_root, case_id)
   %POSTREPAIRSUBLIMATIONCORRELATION Compare complete common daily means.
   manifest_file = fullfile(candidate_root, "eval", "promice", ...
      "manifest.json");
   manifest = jsondecode(fileread(manifest_file));
   cases = manifest.cases;
   selected = string({cases.case_id}) == case_id;
   assert(nnz(selected) == 1, 'PROMICE sign-evidence case is not unique')
   c = cases(selected);
   mar_file = fullfile(candidate_root, "input", "userdata", ...
      string(c.colocation.mar.data_files));
   racmo_file = fullfile(candidate_root, "input", "userdata", ...
      string(c.colocation.racmo.data_files));
   mar = load(mar_file, 'Data');
   racmo = load(racmo_file, 'Data');
   mar_daily = retime(mar.Data(:, "subl"), 'daily', 'mean');
   racmo_daily = retime(racmo.Data(:, "subl"), 'daily', 'mean');
   mar_daily.Properties.VariableNames = "mar_subl";
   racmo_daily.Properties.VariableNames = "racmo_subl";
   paired = synchronize(mar_daily, racmo_daily, 'intersection');
   valid = all(isfinite(paired.Variables), 2);
   assert(nnz(valid) > 365, 'insufficient common sublimation support')
   value = corr(paired.mar_subl(valid), paired.racmo_subl(valid));
end

function rows = repairRows(root, first, second)
   %REPAIRROWS Flatten shared-helper records into durable CSV evidence.
   n = numel(first.records);
   root_kind = repmat(string(root.kind), n, 1);
   artifact_path = strings(n, 1);
   first_status = strings(n, 1);
   second_status = strings(n, 1);
   actions = strings(n, 1);
   changed_variables = strings(n, 1);
   changed_metadata_fields = strings(n, 1);
   hash_before = strings(n, 1);
   hash_after = strings(n, 1);
   second_hash = strings(n, 1);
   unrelated_payload_preserved = false(n, 1);
   for k = 1:n
      artifact_path(k) = root_kind(k) + "/" + string( ...
         erase(first.records(k).filename, string(root.input) + filesep));
      first_status(k) = string(first.records(k).status);
      second_status(k) = string(second.records(k).status);
      actions(k) = strjoin(string(first.records(k).actions), ";");
      changed_variables(k) = strjoin(string( ...
         first.records(k).changed_variables), ";");
      changed_metadata_fields(k) = strjoin(string( ...
         first.records(k).changed_metadata_fields), ";");
      hash_before(k) = string(first.records(k).hash_before);
      hash_after(k) = string(first.records(k).hash_after);
      second_hash(k) = string(second.records(k).hash_after);
      unrelated_payload_preserved(k) = ...
         first.records(k).unrelated_payload_preserved;
   end
   rows = table(root_kind, artifact_path, first_status, second_status, ...
      actions, changed_variables, changed_metadata_fields, hash_before, ...
      hash_after, second_hash, unrelated_payload_preserved);
end

function assertOnlySublimationChanges(report)
   %ASSERTONLYSUBLIMATIONCHANGES Reject any unrelated dry-run transformation.
   allowed_actions = ["flip_racmo_subl_sign", "stamp_racmo_subl_sign", ...
      "restamp_metadata", "sync_artifact_metadata"];
   allowed_metadata = ["racmo_subl_native_sign_convention", ...
      "racmo_subl_sign_convention"];
   for k = 1:numel(report.records)
      assert(all(ismember(string(report.records(k).actions), allowed_actions)), ...
         'bounded repair planned an unrelated action for %s', ...
         report.records(k).filename)
      assert(all(ismember(string( ...
         report.records(k).changed_variables), "subl")), ...
         'bounded repair planned an unrelated value change for %s', ...
         report.records(k).filename)
      assert(all(ismember(string( ...
         report.records(k).changed_metadata_fields), allowed_metadata)), ...
         'bounded repair planned an unrelated metadata change for %s', ...
         report.records(k).filename)
   end
end

function rows = affectedFigureRows(candidate_root)
   %AFFECTEDFIGUREROWS Select exact candidate plots containing RACMO subl.
   filename = fullfile(candidate_root, "preview", "figure_index.csv");
   assert(isfile(filename), 'candidate figure index is missing')
   index = readtable(filename, TextType="string", ...
      VariableNamingRule="preserve");
   has_racmo = hasToken(index.userdata_sources, "racmo2.3p3");
   has_subl = hasToken(index.plotted_variables, "subl");
   rows = index(has_racmo & has_subl, ...
      ["dataset_family", "case_id", "figure_group", "figure_file"]);
   % The candidate was atomically moved after its original figure index was
   % written. Rebuild stable .2.38 promotion-root-relative names instead of
   % preserving the stale absolute root recorded in that historical index.
   rows.figure_file = "candidate/preview/figures/" ...
      + rows.dataset_family + "/" + rows.case_id + "/" ...
      + rows.figure_group + ".png";
end

function tf = hasToken(values, token)
   %HASTOKEN Match one comma-delimited inventory token exactly.
   tf = contains("," + replace(string(values), " ", "") + ",", ...
      "," + token + ",");
end
