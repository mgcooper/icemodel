function runMarProfileBundleProof()
%RUNMARPROFILEBUNDLEPROOF Run bounded MAR profile staging and visual QA.
%
% Run from the repository root with:
%   matlab -batch "icemodel.verification.setup.runMarProfileBundleProof()"
%
% The isolated staging root is deleted on completion. Only the compact proof
% table and three profile figures remain under data/preview/qa.

repo_root = string(icemodel.internal.fullpath());

% Resolve real source caches and allocate a disposable isolated stage.
sumup_dir = icemodel.verification.setup.sumupCacheDir("");
mar_dir = string(getenv("ICEMODEL_MAR_DIR"));
if mar_dir == ""
   mar_dir = "/Volumes/S03/DATA/greenland/mar3p11/RUH2";
end
assert(isfolder(sumup_dir) && isfolder(mar_dir), ...
   'SUMup or MAR source cache is unavailable.')
stage_root = string(tempname);
mkdir(stage_root)
cleanup = onCleanup(@() rmdir(stage_root, 's'));

% These source-backed dates exercise accumulation, upper-ablation/percolation,
% and the KAN_U firn anchor without staging full multi-year forcing windows.
sites = ["egpprofileproof"; "cenprofileproof"; "kanuprofileproof"];
points = [75.62151864466415, -35.973961899732814; ...
   77.15039689818794, -61.062409924023314; ...
   67.00038, -47.02615];
dates = [datetime(2016, 5, 16, TimeZone="UTC"); ...
   datetime(1999, 4, 29, TimeZone="UTC"); ...
   datetime(2016, 4, 26, TimeZone="UTC")];
for n = 1:numel(sites)
   icemodel.verification.setup.importSumup(sumup_dir, ...
      points=points(n, :), case_ids=sites(n), ...
      startdate=dates(n), enddate=dates(n) + days(1) - seconds(1), ...
      output_root=stage_root, build_forcing=true, forcing_sources="mar", ...
      mar_dir=mar_dir, overwrite=true, overwrite_family=false);
end

% Quantitative comparison must resolve the optional sidecar automatically,
% retain exact-date grouping, and produce at least one density match per case.
proof_rows = cell(numel(sites), 7);
output_dir = fullfile(repo_root, 'data', 'preview', 'qa');
proof_dir = fullfile(stage_root, 'proof');
mkdir(proof_dir)
temporary_figures = strings(numel(sites), 1);
final_figures = strings(numel(sites), 1);
relative_figures = strings(numel(sites), 1);
for n = 1:numel(sites)
   result = icemodel.verification.comparecase(sites(n), ...
      evaluation_data_root=fullfile(stage_root, 'eval'), ...
      input_data_root=fullfile(stage_root, 'input'), ...
      dataset_family="sumup", make_plot=false);
   density_rows = result.metrics(result.metrics.variable == "density", :);
   density_status = strjoin(unique(string(density_rows.status)), ',');
   assert(~isempty(density_rows) && any(density_rows.status == "ok"), ...
      ['MAR/SUMup density comparison has no exact-date overlap for %s; ' ...
      'statuses: %s.'], sites(n), density_status)

   temporary_figures(n) = fullfile(proof_dir, ...
      "mar_density_profile_" + sites(n) + ".png");
   final_figures(n) = fullfile(output_dir, ...
       "mar_density_profile_" + sites(n) + ".png");
   % Persist repo-relative evidence paths so the proof ledger is portable.
   relative_figures(n) = fullfile("data", "preview", "qa", ...
       "mar_density_profile_" + sites(n) + ".png");
   f = icemodel.verification.plotcase(sites(n), ...
       evaluation_data_root=fullfile(stage_root, 'eval'), ...
       input_data_root=fullfile(stage_root, 'input'), ...
       dataset_family="sumup", source="compare", visible="off", ...
       output_file=temporary_figures(n));
   close(f)
   proof_rows(n, :) = {sites(n), points(n, 1), points(n, 2), dates(n), ...
       nnz(density_rows.status == "ok"), ...
       min(density_rows.rmse(density_rows.status == "ok")), relative_figures(n)};
end
proof = cell2table(proof_rows, 'VariableNames', { ...
   'case_id', 'lat_wgs84', 'lon_wgs84', 'requested_date', ...
   'ok_density_group_count', 'minimum_density_rmse_kg_m3', 'figure_file'});
temporary_csv = fullfile(proof_dir, 'mar_density_profile_bundle_proof.csv');
final_csv = fullfile(output_dir, 'mar_density_profile_bundle_proof.csv');
writetable(proof, temporary_csv);

% The generic artifact audit must accept the optional model outputs without
% changing primary forcing readiness. The deliberately one-day SUMup target
% tables can retain unrelated native-unit/time findings outside this Bead.
audit = icemodel.verification.auditArtifacts( ...
   evaluation_data_root=fullfile(stage_root, 'eval'), ...
   input_data_root=fullfile(stage_root, 'input'), families="sumup", ...
   report_dir=fullfile(stage_root, 'preview', 'qa'));
profile_artifacts = audit.artifacts( ...
   string({audit.artifacts.kind}) == "model_output");
profile_finding_mask = string({audit.findings.kind}) == "model_output" ...
   & ismember(string({audit.findings.severity}), ["error", "blocker"]);
blocking_findings = audit.findings(profile_finding_mask);
blocking_summary = strjoin(string({blocking_findings.code}) + ": " ...
   + string({blocking_findings.message}), " | ");
assert(numel(profile_artifacts) == numel(sites) ...
   && all([profile_artifacts.exists]) && isempty(blocking_findings), ...
   'Bounded MAR profile model outputs failed artifact QA: %s', blocking_summary)

% Promote only a fully validated proof; a failed run leaves canonical evidence
% untouched and the temporary stage is removed by cleanup.
if ~isfolder(output_dir)
   mkdir(output_dir)
end
for n = 1:numel(sites)
   movefile(temporary_figures(n), final_figures(n), 'f')
end
movefile(temporary_csv, final_csv, 'f')

fprintf('MAR profile proof passed for %d cases; output: %s\n', ...
   height(proof), final_csv)
clear cleanup
end
