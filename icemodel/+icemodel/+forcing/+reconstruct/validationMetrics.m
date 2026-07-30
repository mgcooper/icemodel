function metrics = validationMetrics(truth, filled, gaps, series, channel, kwargs)
   %VALIDATIONMETRICS Grade reconstructed samples against withheld truth.
   %
   %  metrics = icemodel.forcing.reconstruct.validationMetrics( ...
   %     truth, filled, draws.gaps, series, "lwd")
   %
   % Role
   %  Held-out metric computation for the validation protocol (gap-fill
   %  policy). `truth` holds the withheld observed values
   %  and `filled` the method's reconstruction on the same sample axis (NaN
   %  where the method declined); both are restricted to the synthetic-gap
   %  samples by construction. Metrics are reported overall and per gap
   %  stratum (bucket x season), and in-sample statistics never enter here —
   %  callers pass held-out draws only.
   %
   % Name-value
   %  sigma : optional 1-sigma uncertainty per sample (same axis as filled)
   %     for calibration scoring; omit when the method publishes none.
   %  provenance : optional numeric method code per sample. Finite filled
   %     samples must carry a finite non-missing code for complete
   %     provenance accounting.
    %  jump_factor : boundary-jump multiplier over the local median absolute
    %     step (POLICY B6; Section-C default 3).
    %  latitude, longitude : station geometry required for the shortwave
    %     top-of-atmosphere bound.
   %
   % Returns
   %  metrics : struct with fields
   %     overall : one-row table — n, coverage, bias, rmse,
   %        correlation, variability_ratio, within_gap_observed_spread,
   %        bound_violations,
   %        boundary_jump_rate, sigma1_coverage, sigma2_coverage
   %        provenance_accounting, and typical_magnitude (NaN when
   %        undefined or its optional input is absent).
   %     by_stratum : table keyed by bucket and season with the same
   %        measures per stratum.
   %     max_complete_gap_hours : longest held-out gap for which the
   %        candidate supplied every sample (zero when none was complete).
   %
   % See also: icemodel.forcing.reconstruct.syntheticMissingness,
   %  icemodel.forcing.reconstruct.admissionGate,
   %  icemodel.forcing.reconstruct.physicalBounds

   arguments
      truth (:, 1) double
      filled (:, 1) double
      gaps table
      series timetable
      channel (1, 1) string
      kwargs.sigma (:, 1) double = nan(0, 1)
      kwargs.provenance (:, 1) double = nan(0, 1)
       kwargs.jump_factor (1, 1) double {mustBePositive} = ...
          icemodel.forcing.reconstruct.setopts().jump_factor
       kwargs.latitude (1, 1) double = NaN
       kwargs.longitude (1, 1) double = NaN
   end
   if numel(truth) ~= numel(filled)
      error('icemodel:reconstruct:validationMetrics:sizeMismatch', ...
         'truth and filled must share one sample axis');
   end
   if ~isempty(kwargs.provenance) ...
         && numel(kwargs.provenance) ~= numel(filled)
      error('icemodel:reconstruct:validationMetrics:provenanceSizeMismatch', ...
         'provenance and filled must share one sample axis');
   end
   if ~issorted(gaps.start_time)
      % The sequential sample-to-gap cursor below assumes time order, which
      % syntheticMissingness guarantees; refuse anything else loudly.
      error('icemodel:reconstruct:validationMetrics:unsortedGaps', ...
         'gaps must be sorted by start_time');
   end

   times = series.Properties.RowTimes;
   x = series.(channel);
    % The local step scale contextualizes boundary jumps: POLICY B6 keys
   % it to the station AND season, through the shared helper so the engine
   % tiers reject exactly what these metrics would score as violations.
   season_scale = icemodel.forcing.reconstruct.stepScale(times, x);

   % Map each truth/filled sample to its gap row so strata aggregate cleanly.
    sample_gap = zeros(numel(truth), 1);
    sample_index = zeros(numel(truth), 1);
   sample_cursor = 0;
   jump_violation = false(height(gaps), 1);
   for g = 1:height(gaps)
       idx = find(times >= gaps.start_time(g) & times <= gaps.end_time(g));
       sample_gap(sample_cursor + 1:sample_cursor + numel(idx)) = g;
       sample_index(sample_cursor + 1:sample_cursor + numel(idx)) = idx;
      sample_cursor = sample_cursor + numel(idx);
      % Boundary jump: the step between the last native sample before the
      % gap (and first after) and the adjacent filled sample.
      first_fill = filled(sample_cursor - numel(idx) + 1);
      last_fill = filled(sample_cursor);
      before = idx(1) - 1;
      after = idx(end) + 1;
      jumps = nan(2, 1);
      if before >= 1 && isfinite(x(before)) && isfinite(first_fill)
         jumps(1) = abs(first_fill - x(before));
      end
      if after <= numel(x) && isfinite(x(after)) && isfinite(last_fill)
         jumps(2) = abs(last_fill - x(after));
      end
       measurable = isfinite(jumps);
       jump_violation(g) = ~any(measurable) ...
          || any(jumps(measurable) > kwargs.jump_factor * ...
          season_scale.(char(gaps.season(g))));
   end
    if sample_cursor ~= numel(truth)
      error('icemodel:reconstruct:validationMetrics:gapSampleMismatch', ...
         'gap table spans %d samples but truth has %d', ...
         sample_cursor, numel(truth));
    end

   % Duration authorization needs a complete held-out gap, not merely a
   % finite fragment somewhere inside a longer draw. Partial coverage still
   % contributes to the skill metrics and can admit a donor when other draws
   % of the same stratum are complete.
   complete_gap = false(height(gaps), 1);
   for g = 1:height(gaps)
      complete_gap(g) = all(isfinite(filled(sample_gap == g)));
   end
   max_complete_gap_hours = max([gaps.duration_hours(complete_gap); 0]);

    % One shared validator makes metrics and production reject the same
    % scalar and relational violations.
    sample_times = times(sample_index);
    if channel == "swu"
       physical_valid = icemodel.forcing.reconstruct.physicalValidity( ...
          channel, filled, sample_times, swd=series.swd(sample_index));
    else
       physical_valid = icemodel.forcing.reconstruct.physicalValidity( ...
          channel, filled, sample_times, latitude=kwargs.latitude, ...
          longitude=kwargs.longitude, interval=median(diff(times)));
    end

     overall = scoreSubset(truth, filled, kwargs.sigma, kwargs.provenance, ...
       physical_valid, jump_violation, sample_gap, ...
       true(numel(truth), 1), ...
       true(height(gaps), 1));

   % Per-stratum rows over the bucket x season combinations present.
   [strata, ~, stratum_of_gap] = unique(gaps(:, {'bucket', 'season'}), 'rows');
   stratum_rows = cell(height(strata), 1);
   for s = 1:height(strata)
      gap_in = stratum_of_gap == s;
      sample_in = ismember(sample_gap, find(gap_in));
        stratum_rows{s} = scoreSubset(truth, filled, kwargs.sigma, ...
           kwargs.provenance, physical_valid, ...
           jump_violation, sample_gap, sample_in, gap_in);
   end
   by_stratum = [strata, vertcat(stratum_rows{:})];

   metrics = struct('overall', overall, 'by_stratum', by_stratum, ...
      'max_complete_gap_hours', max_complete_gap_hours);
end

function row = scoreSubset(truth, filled, sigma, provenance, physical_valid, ...
      jump_violation, sample_gap, sample_in, gap_in)
   %SCORESUBSET Compute one metric row over a sample/gap subset.
   t = truth(sample_in);
   f = filled(sample_in);
   have = isfinite(f);
   err = f(have & isfinite(t)) - t(have & isfinite(t));
   paired = have & isfinite(t);
   observed = t(paired);
   predicted = f(paired);
   typical_magnitude = mean(abs(observed), 'omitnan');

   % RMSE and bias alone reward regression toward the mean. Correlation uses
   % all paired samples, while variability removes each held-out gap's mean
   % first so between-gap seasonal offsets cannot hide weather compression.
   correlation = NaN;
   variability_ratio = NaN;
   within_gap_observed_spread = NaN;
   if numel(observed) >= 2
      predicted_spread = std(predicted);
      if predicted_spread > 0
         C = corrcoef(observed, predicted);
         correlation = C(1, 2);
      end
   end
   gap_id = sample_gap(sample_in);
   observed_anomaly = nan(size(t));
   predicted_anomaly = nan(size(f));
   complete_spread_evidence = true;
   for g = reshape(find(gap_in), 1, [])
      in_gap = paired & gap_id == g;
      if nnz(in_gap) < 2
         complete_spread_evidence = false;
         continue
      end
      observed_anomaly(in_gap) = t(in_gap) - mean(t(in_gap));
      predicted_anomaly(in_gap) = f(in_gap) - mean(f(in_gap));
   end
   anomaly_pair = isfinite(observed_anomaly) & isfinite(predicted_anomaly);
   if complete_spread_evidence && nnz(anomaly_pair) >= 2
      within_gap_observed_spread = std(observed_anomaly(anomaly_pair));
      if within_gap_observed_spread > 0
         variability_ratio = std(predicted_anomaly(anomaly_pair)) ...
            / within_gap_observed_spread;
      end
   end

   % Uncertainty calibration: empirical coverage of nominal 1/2-sigma
   % intervals; NaN when the method publishes no uncertainty.
   sigma1 = NaN;
   sigma2 = NaN;
   if ~isempty(sigma)
      s = sigma(sample_in);
      ok = have & isfinite(t) & isfinite(s) & s > 0;
      if any(ok)
         z = abs(filled(sample_in) - truth(sample_in)) ./ sigma(sample_in);
         z = z(ok);
         sigma1 = mean(z <= 1);
         sigma2 = mean(z <= 2);
      end
   end

   % Provenance is graded only for values a candidate actually supplied;
   % declined samples do not need a method code.
   provenance_accounting = NaN;
   if ~isempty(provenance) && any(have)
      p = provenance(sample_in);
      provenance_accounting = mean(isfinite(p(have)) & p(have) ~= 255);
   end

     valid = physical_valid(sample_in);
     row = table(nnz(sample_in), mean(have), mean(err), sqrt(mean(err.^2)), ...
        correlation, variability_ratio, within_gap_observed_spread, ...
        typical_magnitude, ...
        nnz(have & ~valid), mean(jump_violation(gap_in)), sigma1, sigma2, ...
        provenance_accounting, ...
       'VariableNames', {'n', 'coverage', 'bias', 'rmse', ...
       'correlation', 'variability_ratio', 'within_gap_observed_spread', ...
       'typical_magnitude', ...
       'bound_violations', 'boundary_jump_rate', 'sigma1_coverage', ...
      'sigma2_coverage', 'provenance_accounting'});
end
