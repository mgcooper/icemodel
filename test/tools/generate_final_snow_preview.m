%GENERATE_FINAL_SNOW_PREVIEW Build canonical seasonal QA, figures, and readiness.
% Reads promoted data/input and data/eval artifacts, then writes only the
% canonical data/preview products consumed by the checked Quarto report.

driver_root = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(driver_root));
data_root = fullfile(repo_root, "data");
eval_root = fullfile(data_root, "eval");
input_root = fullfile(data_root, "input");
preview_root = fullfile(data_root, "preview");
figure_root = fullfile(preview_root, "figures");
qa_root = fullfile(preview_root, "qa");
evidence_root = fullfile(repo_root, "sandbox", "verification", ...
   "icemodel-hfc.2.38", "evidence");
complete_file = fullfile(evidence_root, "FINAL_PREVIEW_COMPLETE.json");
failed_file = fullfile(evidence_root, "FINAL_PREVIEW_FAILED.txt");
promotion_file = fullfile(evidence_root, "COMPLETE.json");

% Final generation is invalid until the manifest-backed data promotion has
% converged and published its completion marker.
assert(isfile(promotion_file), ...
   "icemodel:hfc238:promotionIncomplete", ...
   "Run and verify the seasonal artifact promotion before final preview generation.");
if ~isfolder(evidence_root)
   mkdir(evidence_root)
end
if isfile(complete_file)
   delete(complete_file)
end
if isfile(failed_file)
   delete(failed_file)
end

% Resolve the repository package and external plotting dependencies once.
addpath(fullfile(repo_root, "icemodel"));
icemodel.dependencies();

try
   started = datetime("now", TimeZone="UTC");
   families = ["promice", "esm_snowmip", "laugh_tests"];

   % The promoted manifests must be self-contained release products; raw source
   % paths may be external, but no selected manifest may point to a temp stage.
   for family = families
      manifest_file = fullfile(eval_root, family, "manifest.json");
      assert(isfile(manifest_file), ...
         "icemodel:hfc238:manifestMissing", ...
         "Missing promoted manifest: %s", manifest_file);
      assert(~contains(fileread(manifest_file), "/private/tmp"), ...
         "icemodel:hfc238:temporaryManifestPath", ...
         "Promoted manifest still depends on /private/tmp: %s", manifest_file);
   end

   % First-class QA runs before rendering so actionable artifact defects do not
   % waste a complete 62-case figure pass.
   if ~isfolder(qa_root)
      mkdir(qa_root)
   end
   audit = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      families=families, report_dir=qa_root);
   % NUK_K is explicitly RACMO-unavailable after the native ice-mask audit,
   % reducing the accepted seasonal inventory by one referenced artifact.
   assert(audit.summary.passed && audit.summary.error_count == 0 ...
      && audit.summary.blocker_count == 0 ...
      && audit.summary.case_count == 62 ...
      && audit.summary.artifact_count == 347, ...
      "icemodel:hfc238:artifactAuditFailed", ...
      "Seasonal artifact audit failed or its 62-case/347-artifact scope drifted.");

   % Recompute the handoff from the exact promoted manifest rather than copying
   % the earlier candidate tables.
   readiness = icemodel.verification.setup.writePromiceSnowModelReadyYears( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      output_dir=qa_root);
   check = readiness.check;
   assert(check.all_checks_passed ...
      && check.candidate_site_year_count == 550 ...
      && check.ready_site_year_count == 115 ...
      && check.strict_site_year_count == 107 ...
      && check.ready_site_count == 22, ...
      "icemodel:hfc238:readinessMismatch", ...
      "Promoted PROMICE readiness counts differ from the validated handoff.");

   % Every accepted figure is reproducible, so remove only the three owned
   % family trees and rebuild them with the current renderer and source colors.
   for family = families
      family_root = fullfile(figure_root, family);
      if isfolder(family_root)
         rmdir(family_root, 's')
      end
   end
   if ~isfolder(figure_root)
      mkdir(figure_root)
   end
summary = icemodel.verification.plotVerificationArtifacts( ...
   dataset_family=families, output_root=data_root, ...
   figure_root=figure_root, save_figs=true, overwrite=true, visible=false);

% Firn previews share the canonical figure root but are independently owned.
% Count only the three seasonal family trees rebuilt by this driver so a later
% current-code refresh cannot mistake accepted firn PNGs for seasonal drift.
pngs = dir(fullfile(figure_root, families(1), "**", "*.png"));
for family_idx = 2:numel(families)
   pngs = [pngs; dir(fullfile( ...
      figure_root, families(family_idx), "**", "*.png"))]; %#ok<AGROW>
end
assert(istable(summary) && height(summary) == numel(pngs) ...
   && height(summary) == 678, ...
      "icemodel:hfc238:figureInventoryMismatch", ...
      "Expected the validated 678-row/PNG final figure inventory.");
   writetable(summary, fullfile(preview_root, "figure_index.csv"));
   save(fullfile(qa_root, "final_snow_preview_summary.mat"), ...
      "summary", "check");

   % Publish a small completion record; Quarto generation and its independent
   % link/freshness checker run afterward from repository-owned Python tools.
   completed = datetime("now", TimeZone="UTC");
   evidence = struct( ...
      'started_at', string(started), 'completed_at', string(completed), ...
      'elapsed_seconds', seconds(completed - started), ...
      'figure_count', numel(pngs), 'summary_rows', height(summary), ...
      'audit_error_count', audit.summary.error_count, ...
      'audit_blocker_count', audit.summary.blocker_count, ...
      'ready_site_year_count', check.ready_site_year_count, ...
      'strict_site_year_count', check.strict_site_year_count);
   writelines(string(jsonencode(evidence, PrettyPrint=true)), complete_file);
catch err
   writelines(string(getReport(err, "extended", "hyperlinks", "off")), ...
      failed_file);
   rethrow(err)
end
