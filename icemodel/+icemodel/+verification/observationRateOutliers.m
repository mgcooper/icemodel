function rates = observationRateOutliers(summary, policy)
   %OBSERVATIONRATEOUTLIERS Flag site-years whose observed ablation rate is far
   %below the same station's own family of scored site-years.
   %
   % An observation record can be systematically compressed without tripping
   % any readiness gate, because the gates only reject FLAGGED transitions and
   % unresolved steps. An unflagged sensor or datum problem passes admission
   % and appears in the comparison as a model error. Comparing each station
   % against its own distribution catches that without penalising genuinely
   % low-melt sites.
   %
   % The flag is a caveat, never an exclusion. Dropping the flagged site-years
   % would raise apparent model skill by removing observations, so the decision
   % to exclude stays with the reader.

   completed = string(summary.status) == "completed";
   case_id = string(summary.case_id(completed));
   site_id = string(summary.site_id(completed));
   year = double(summary.year(completed));
   observed_mwe = double(summary.observation_intact_mwe(completed));

   % The window length varies between site-years, so a rate is the only
   % comparable quantity across a station's record.
   window_days = days(summary.window_end(completed) ...
      - summary.window_start(completed));
   observed_rate = observed_mwe ./ window_days;
   observed_rate(~isfinite(window_days) | window_days <= 0) = NaN;

   % Each station is scored against its own median, so a station that simply
   % melts less than its neighbours is not flagged for that alone. The family
   % counts completed site-years that produced a finite rate: a case with no
   % usable window contributes no rate and so cannot enlarge the family.
   station_median_rate = NaN(size(observed_rate));
   station_scored_year_count = zeros(size(observed_rate));
   stations = unique(site_id, 'stable');
   for k = 1:numel(stations)
      member = site_id == stations(k);
      family = observed_rate(member & isfinite(observed_rate));
      station_scored_year_count(member) = numel(family);
      if ~isempty(family)
         station_median_rate(member) = median(family);
      end
   end

   % Dividing by a negative median would flip the comparison and flag the
   % station's HIGHEST-ablation year. A station that nets accumulation over
   % its admitted years has no meaningful ablation-rate family.
   rate_ratio = NaN(size(observed_rate));
   usable_median = station_median_rate > 0;
   rate_ratio(usable_median) = observed_rate(usable_median) ...
      ./ station_median_rate(usable_median);
   % A family needs enough members for its median to mean anything.
   has_family = station_scored_year_count ...
      >= policy.observation_rate_outlier_min_years;
   flagged = has_family & isfinite(rate_ratio) ...
      & rate_ratio < policy.observation_rate_outlier_ratio;

   rates = table(case_id, site_id, year, observed_mwe, window_days, ...
      observed_rate, station_median_rate, station_scored_year_count, ...
      rate_ratio, ...
      flagged, 'VariableNames', {'case_id', 'site_id', 'year', ...
      'observation_intact_mwe', 'window_days', ...
      'observed_rate_mwe_per_day', 'station_median_rate_mwe_per_day', ...
      'station_scored_year_count', 'rate_ratio_to_station_median', ...
      'rate_outlier_flag'});
   rates = sortrows(rates, {'site_id', 'year'});
end
