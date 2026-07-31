function skill = commonSupportSkill(truth, candidate, baseline)
   %COMMONSUPPORTSKILL Compare candidate and baseline on identical samples.
   %
   %  skill = icemodel.forcing.reconstruct.commonSupportSkill( ...
   %     truth, candidate, baseline)
   %
   % Candidate coverage is graded separately. This helper prevents a
   % partially supported method from claiming skill by selecting samples
   % that are easier than the samples used to score its baseline.

   arguments
      truth (:, 1) double
      candidate (:, 1) double
      baseline (:, 1) double
   end

   if numel(candidate) ~= numel(truth) || numel(baseline) ~= numel(truth)
      error('icemodel:reconstruct:commonSupportSkill:sizeMismatch', ...
         'truth, candidate, and baseline must share one sample axis');
   end

   common = isfinite(truth) & isfinite(candidate) & isfinite(baseline);
   if ~any(common)
      skill = struct('n', 0, 'candidate_rmse', NaN, ...
         'baseline_rmse', NaN, 'rmse_ratio', NaN, ...
         'fractional_improvement', NaN);
      return
   end

   candidate_rmse = sqrt(mean((candidate(common) - truth(common)).^2));
   baseline_rmse = sqrt(mean((baseline(common) - truth(common)).^2));
   if baseline_rmse == 0
      % Exact ties have zero improvement; a nonzero candidate is infinitely
      % worse than an exact baseline.
      if candidate_rmse == 0
         rmse_ratio = 1;
      else
         rmse_ratio = Inf;
      end
   else
      rmse_ratio = candidate_rmse / baseline_rmse;
   end
   skill = struct('n', nnz(common), 'candidate_rmse', candidate_rmse, ...
      'baseline_rmse', baseline_rmse, 'rmse_ratio', rmse_ratio, ...
      'fractional_improvement', 1 - rmse_ratio);
end
