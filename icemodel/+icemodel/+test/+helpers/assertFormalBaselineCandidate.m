function assertFormalBaselineCandidate(kind, candidate, cases, selector)
   %ASSERTFORMALBASELINECANDIDATE Reject incomplete state before publication.
   %
   %  icemodel.test.helpers.assertFormalBaselineCandidate( ...
   %     "regression", candidate, cases, "rolling")

   arguments
      kind (1, 1) string {mustBeMember(kind, ["regression", "perf"])}
      candidate table
      cases table
      selector (1, 1) string
   end

   % Candidate membership must be an exact one-to-one image of the selected
   % formal case matrix before any predecessor is archived.
   if isempty(candidate) || ~ismember('case_id', ...
         candidate.Properties.VariableNames)
      error('icemodel:test:baselineCandidateCasesInvalid', ...
         'Formal baseline candidate has no case identity rows.')
   end
   expected = icemodel.test.helpers.normalizeFormalCaseId(string(cases.case_id));
   actual = icemodel.test.helpers.normalizeFormalCaseId(string(candidate.case_id));
   if height(candidate) ~= height(cases) ...
         || numel(unique(actual)) ~= numel(actual) ...
         || ~isequal(sort(actual), sort(expected))
      error('icemodel:test:baselineCandidateCasesInvalid', ...
         'Formal baseline candidate does not exactly match selected cases.')
   end
   icemodel.test.helpers.assertFormalBaselineForcing(candidate, selector);

   % Numerical candidates require finite generated metrics; performance
   % candidates additionally require valid samples and positive timing state.
   switch kind
      case "regression"
         % Closed metadata allowlist. baseline_tag is a persisted build
         % identity, not a numerical regression metric; every column not
         % listed here is validated as a metric.
         metadata = ["case_id", "tier", "baseline_type", "baseline_tag", ...
            "smbmodel", "sitename", "forcings", "simyear", "solver", ...
            "last_updated_utc"];
         metrics = setdiff(string(candidate.Properties.VariableNames), ...
            metadata, 'stable');
         invalid = isempty(metrics);
         for metric = metrics
            values = candidate.(metric);
            invalid = invalid || ~isnumeric(values) ...
               || any(~isfinite(values), 'all');
         end
      case "perf"
         required = ["valid", "passed_perf", "median_wall_s", ...
            "mean_wall_s", "min_wall_s", "max_wall_s", "n_runs", ...
            "n_warmups", "tol_perf"];
         names = string(candidate.Properties.VariableNames);
         invalid = ~all(ismember(required, names));
         if ~invalid
            positive = [candidate.median_wall_s, candidate.mean_wall_s, ...
               candidate.min_wall_s, candidate.max_wall_s, ...
               candidate.n_runs, candidate.tol_perf];
            invalid = any(~logical(candidate.valid)) ...
               || any(~logical(candidate.passed_perf)) ...
               || any(~isfinite(positive) | positive <= 0, 'all') ...
               || any(~isfinite(candidate.n_warmups) ...
               | candidate.n_warmups < 0);
         end
   end
   if invalid
      error('icemodel:test:baselineCandidateMetricsInvalid', ...
         'Formal %s baseline candidate contains incomplete or invalid evidence.', ...
         kind)
   end
end
