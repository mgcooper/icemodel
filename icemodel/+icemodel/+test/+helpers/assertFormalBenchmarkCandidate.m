function assertFormalBenchmarkCandidate(candidate)
   %ASSERTFORMALBENCHMARKCANDIDATE Reject invalid component timing evidence.
   required = ["Valid", "SampleSize", "Mean"];
   names = string(candidate.Properties.VariableNames);
   invalid = isempty(candidate) || ~all(ismember(required, names));
   if ~invalid
      invalid = any(~logical(candidate.Valid)) ...
         || any(~isfinite(candidate.SampleSize) | candidate.SampleSize <= 0) ...
         || any(~isfinite(candidate.Mean) | candidate.Mean <= 0);
   end
   if invalid
      error('icemodel:test:benchmarkBaselineCandidateInvalid', ...
         'Benchmark baseline candidate contains incomplete or invalid evidence.')
   end
end
