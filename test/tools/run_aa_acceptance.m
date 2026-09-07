function report = run_aa_acceptance(artifacts_a, artifacts_b)
   %RUN_AA_ACCEPTANCE Check that the formal timing protocol reproduces itself.
   %
   %  report = run_aa_acceptance()
   %  report = run_aa_acceptance(artifacts_a, artifacts_b)
   %
   % The A/A test runs the same code twice under the formal process-isolated
   % protocol and requires the two runs to reproduce each other. It validates
   % the measurement system, not the code. Run it before trusting an A/B
   % comparison.
   %
   % With no arguments, run both measurement passes back to back and compare
   % them. With arguments, compare already-saved artifacts and run nothing
   % (use this form when the two runs were made by hand). Each argument is
   % one artifact path or a string array of them: a run with the default
   % smbmodel="all" writes one artifact per model, so pass both files from
   % run A and both from run B. Pairing is by sorted filename, which aligns
   % the per-model artifacts.
   %
   % Pass criteria:
   %  - each side has the same artifact count and no duplicate paths
   %  - the two sides share no paths
   %  - each side has one meta.run_name, and the two names differ
   %  - every artifact has the same nonempty meta.hostname
   %  - every artifact has the same MATLAB version and input data root
   %  - every artifact has the same nonempty meta.git_revision
   %  - tier, simyear, n_runs, n_warmups, and tol_perf match
   %  - the case ID and forcing product match for every measured case
   %  - every run artifact reports meta.isolation = "process" (the only
   %    protocol this gate certifies) and meta.ambient_stable = true
   %  - every joined case reports valid = true in both runs
   %  - every per-case median ratio B/A lies inside [1/AA_BAND, AA_BAND]
   %
   % Output report fields: passed, ambient_stable, aa_band, and a per-case
   % table with both medians and their ratio.

   if nargin == 0
      % Two full formal runs, back to back, in this session's process. Each
      % measured case still runs in its own fresh subprocess, so the two
      % passes share nothing but the host.
      results_a = run_perf_suite(isolation="process", ...
         build_report=false, include_benchmarks=false);
      results_b = run_perf_suite(isolation="process", ...
         build_report=false, include_benchmarks=false);
      artifacts_a = string(results_a.artifact_file);
      artifacts_b = string(results_b.artifact_file);
   end

   % Sort by filename so parent directories do not change the model pairing.
   artifacts_a = reshape(string(artifacts_a), [], 1);
   artifacts_b = reshape(string(artifacts_b), [], 1);
   [~, names_a, ext_a] = fileparts(artifacts_a);
   [~, names_b, ext_b] = fileparts(artifacts_b);
   [~, order_a] = sort(names_a + ext_a);
   [~, order_b] = sort(names_b + ext_b);
   artifacts_a = artifacts_a(order_a);
   artifacts_b = artifacts_b(order_b);
   assert(numel(artifacts_a) == numel(artifacts_b), ...
      'icemodel:test:aaAcceptance:artifactCountMismatch', ...
      'run A and run B must supply the same number of artifacts')
   assert(numel(unique(artifacts_a)) == numel(artifacts_a) ...
      && numel(unique(artifacts_b)) == numel(artifacts_b), ...
      'icemodel:test:aaAcceptance:duplicateArtifact', ...
      'each run must supply each artifact once')

   % One run compared against itself reproduces trivially, so the two sides
   % must not share a file.
   assert(~any(ismember(artifacts_a, artifacts_b)), ...
      'icemodel:test:aaAcceptance:reusedFile', ...
      'run_aa_acceptance needs two independent runs, got the same file')

   % Load every artifact once; the side checks below need all of them.
   A_all = arrayfun(@load, artifacts_a, 'UniformOutput', false);
   B_all = arrayfun(@load, artifacts_b, 'UniformOutput', false);

   % Each side must be one run: a side that mixes runs could pair
   % artifacts whose differences cancel in the ratios. meta.run_name is
   % the date-qualified identity (yyyymmdd-HHMMSS); meta.run_id holds
   % the time only, which repeats across days.
   run_names_a = cellfun(@(s) string(s.meta.run_name), A_all);
   run_names_b = cellfun(@(s) string(s.meta.run_name), B_all);
   assert(all(~ismissing(run_names_a)) ...
      && all(strlength(strip(run_names_a)) > 0) ...
      && all(~ismissing(run_names_b)) ...
      && all(strlength(strip(run_names_b)) > 0), ...
      'icemodel:test:aaAcceptance:unknownRunName', ...
      'run_aa_acceptance needs a recorded run_name in every artifact')
   assert(isscalar(unique(run_names_a)) ...
      && isscalar(unique(run_names_b)), ...
      'icemodel:test:aaAcceptance:mixedSide', ...
      'run_aa_acceptance needs one run per side, got mixed run_names')
   assert(run_names_a(1) ~= run_names_b(1), ...
      'icemodel:test:aaAcceptance:reusedRun', ...
      'run_aa_acceptance needs two independent runs, got one run_name')

   % Timings from different machines or MATLAB versions are not
   % comparable (see icemodel.test.helpers.perfBaselineCompatibility).
   % meta.hostname names the machine; meta.host only names the platform.
   hosts = [cellfun(@(s) string(s.meta.hostname), A_all); ...
      cellfun(@(s) string(s.meta.hostname), B_all)];
   versions = [cellfun(@(s) string(s.meta.matlab_version), A_all); ...
      cellfun(@(s) string(s.meta.matlab_version), B_all)];
   assert(all(~ismissing(hosts)) && all(strlength(strip(hosts)) > 0), ...
      'icemodel:test:aaAcceptance:unknownHost', ...
      'run_aa_acceptance needs a recorded hostname in every artifact')
   assert(all(~ismissing(versions)) && all(strlength(strip(versions)) > 0), ...
      'icemodel:test:aaAcceptance:unknownMatlabVersion', ...
      'run_aa_acceptance needs a recorded MATLAB version in every artifact')
   assert(isscalar(unique(hosts)) && isscalar(unique(versions)), ...
      'icemodel:test:aaAcceptance:environmentMismatch', ...
      'run_aa_acceptance needs one machine and MATLAB version for both runs')

   % Two runs of different code certify an A/B change, not reproducibility.
   % The revision contains the Git description and a content digest for a
   % dirty tree. Two unchanged back-to-back runs therefore compare equal.
   revisions = [cellfun(@(s) string(s.meta.git_revision), A_all); ...
      cellfun(@(s) string(s.meta.git_revision), B_all)];
   assert(all(~ismissing(revisions)) ...
      && all(strlength(strip(revisions)) > 0), ...
      'icemodel:test:aaAcceptance:unknownRevision', ...
      'run_aa_acceptance needs a recorded git revision in every artifact')
   assert(isscalar(unique(revisions)), ...
      'icemodel:test:aaAcceptance:revisionMismatch', ...
      'run_aa_acceptance needs one source revision for both runs')

   % Medians from different measurement procedures are not comparable:
   % a three-sample run must not certify a thirty-sample protocol. The
   % numeric fields compare exactly, not through a formatted string.
   all_runs = [A_all; B_all];
   tiers = cellfun(@(s) string(s.meta.tier), all_runs);
   simyears = cellfun(@(s) s.meta.simyear, all_runs);
   sample_counts = cellfun(@(s) s.meta.n_runs, all_runs);
   warmup_counts = cellfun(@(s) s.meta.n_warmups, all_runs);
   tolerances = cellfun(@(s) s.meta.tol_perf, all_runs);
   assert(isscalar(unique(tiers)) && isscalar(unique(simyears)) ...
      && isscalar(unique(sample_counts)) ...
      && isscalar(unique(warmup_counts)) ...
      && isscalar(unique(tolerances)), ...
      'icemodel:test:aaAcceptance:procedureMismatch', ...
      'run_aa_acceptance needs one measurement procedure for both runs')

   % Runs that measured different input trees are not comparable even
   % when every case_id matches.
   data_roots = cellfun(@(s) string(s.meta.data_root), all_runs);
   assert(all(~ismissing(data_roots)) ...
      && all(strlength(strip(data_roots)) > 0), ...
      'icemodel:test:aaAcceptance:unknownDataRoot', ...
      'run_aa_acceptance needs a recorded input data root in every artifact')
   assert(isscalar(unique(data_roots)), ...
      'icemodel:test:aaAcceptance:inputMismatch', ...
      'run_aa_acceptance needs one input data root for both runs')

   % The band matches the runs' own two-sided noise budget: artifacts
   % made under a nondefault tol_perf certify against that tolerance,
   % not the current perfMeasurementPolicy default.
   AA_BAND = 1 + A_all{1}.meta.tol_perf;

   % Compare each artifact pair; collect the per-pair tables and stack once.
   ambient_stable = true;
   pair_cases = cell(numel(artifacts_a), 1);
   for k = 1:numel(artifacts_a)
      A = A_all{k};
      B = B_all{k};

      % This gate certifies the process-isolated protocol only; a
      % session-isolated artifact would certify a different protocol.
      assert(A.meta.isolation == "process" ...
         && B.meta.isolation == "process", ...
         'icemodel:test:aaAcceptance:notProcessIsolated', ...
         'run_aa_acceptance compares process-isolated artifacts only')

      % An anchor-invalid run cannot certify anything.
      ambient_stable = ambient_stable ...
         && A.meta.ambient_stable && B.meta.ambient_stable;

      % Match the case ID and forcing product. A case ID does not contain
      % the forcing product.
      a = A.case_summary;
      b = B.case_summary;
      assert(~isempty(a) && ~isempty(b), ...
         'icemodel:test:aaAcceptance:emptyCaseSummary', ...
         'every A/A artifact must contain at least one measured case')
      identity_a = [string(a.case_id), string(a.forcings)];
      identity_b = [string(b.case_id), string(b.forcings)];
      [~, ia, ib] = intersect(identity_a, identity_b, 'rows', 'stable');
      assert(numel(ia) == height(a) && numel(ib) == height(b), ...
         'icemodel:test:aaAcceptance:caseIdentityMismatch', ...
         'run A and run B measured different case workloads')
      ratio = b.median_wall_s(ib) ./ a.median_wall_s(ia);

      % An invalid measurement saves a median anyway; a ratio of two
      % noisy medians must not certify the protocol.
      measurement_valid = a.valid(ia) & b.valid(ib);
      pair_cases{k} = table(a.case_id(ia), a.median_wall_s(ia), ...
         b.median_wall_s(ib), ratio, measurement_valid, 'VariableNames', ...
         {'case_id', 'median_a_s', 'median_b_s', 'ratio_b_over_a', ...
         'measurement_valid'});
   end
   cases = vertcat(pair_cases{:});

   in_band = cases.ratio_b_over_a >= 1 / AA_BAND ...
      & cases.ratio_b_over_a <= AA_BAND;
   passed = ambient_stable && all(in_band) ...
      && all(cases.measurement_valid);

   report = struct('passed', passed, 'aa_band', AA_BAND, ...
      'ambient_stable', ambient_stable, 'cases', cases);

   % Print the verdict and the table so a GUI or -batch caller sees both.
   if passed
      fprintf('A/A PASSED: all %d case ratios inside [%.3f, %.3f]\n', ...
         height(cases), 1 / AA_BAND, AA_BAND);
   else
      fprintf('A/A FAILED: ambient_stable=%d in_band=%d/%d\n', ...
         ambient_stable, sum(in_band), height(cases));
   end
   disp(cases)
end
