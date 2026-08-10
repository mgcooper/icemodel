function [first, last] = evaluationSeason(y, policy)
   %EVALUATIONSEASON Summertime display and evaluation bounds for one year.
   %
   %  [first, last] = icemodel.verification.helpers.evaluationSeason(y, policy)
   %
   % Readiness admits a site-year against this season and the runner evaluates
   % it against the same season.
   %
   % Inputs
   %  y      - calendar year
   %  policy - PROMICE ablation policy carrying the season month/day bounds
   %
   % Outputs
   %  first - season start, UTC
   %  last  - season end, UTC

   start_md = policy.evaluation_season_start_month_day;
   end_md = policy.evaluation_season_end_month_day;
   first = datetime(y, start_md(1), start_md(2), 'TimeZone', 'UTC');
   last = datetime(y, end_md(1), end_md(2), 'TimeZone', 'UTC');
end
