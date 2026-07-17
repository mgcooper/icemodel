function report = reportPromiceCoverage(coverage, requested, leg_windows)
   %REPORTPROMICECOVERAGE Print requested-vs-actual source coverage.
   %
   %  report = icemodel.verification.setup.reportPromiceCoverage( ...
   %     coverage, requested, leg_windows)
   %
   % Renders a per-source coverage table comparing the requested study window
   % against what each forcing source actually delivers on disk, listing the
   % missing years explicitly. Honest reporting is the point: the user is told
   % exactly where the requested window is NOT covered.
   %
   % Inputs
   %  coverage     struct from icemodel.verification.setup.promiceSourceCoverage
   %  requested    [year_start year_end] requested study window (years)
   %  leg_windows  struct keyed by model, each with .start/.end datetime (the
   %               actual per-leg window the driver will stage), or "skipped".
   %
   % Outputs
   %  report  string : the rendered table (also printed to the command window)
   %
   % See also: icemodel.verification.setup.promiceSourceCoverage,
   %  icemodel.verification.setup.importPromiceSites

   arguments
      coverage (1, 1) struct
      requested (1, 2) double
      leg_windows (1, 1) struct = struct()
   end

   y0 = requested(1);
   y1 = requested(2);
   requested_years = y0:y1;

   lines = strings(0, 1);
   lines(end+1) = sprintf('PROMICE co-location source coverage (requested %d-%d)', y0, y1);
   lines(end+1) = string(repmat('-', 1, 78));
   lines(end+1) = sprintf('%-8s %-14s %-22s %s', ...
      'source', 'on-disk', 'staged window', 'missing vs requested');
   lines(end+1) = string(repmat('-', 1, 78));

   % PROMICE anchor is always full-record; report it as covering the request.
   lines(end+1) = sprintf('%-8s %-14s %-22s %s', ...
      'promice', 'full record', stagedTag(leg_windows, 'promice', y0, y1), 'none');

   sources = icemodel.verification.namelists.rcmsources();
   for m = sources
      if ~isfield(coverage, m)
         continue
      end
      cov = coverage.(m);
      if isempty(cov.years)
         disk = "absent";
         missing = sprintf('ALL (%s)', cov.reason);
      else
         disk = sprintf('%d-%d', cov.year_min, cov.year_max);
         miss = setdiff(requested_years, cov.years);
         missing = compactYears(miss);
      end
      lines(end+1) = sprintf('%-8s %-14s %-22s %s', ...
         m, disk, stagedTag(leg_windows, char(m), y0, y1), missing); %#ok<AGROW>
   end

   lines(end+1) = string(repmat('-', 1, 78));
   report = strjoin(lines, newline);
   fprintf('%s\n', report);
end

%% Local functions
function tag = stagedTag(leg_windows, model, y0, y1)
   if ~isfield(leg_windows, model)
      tag = sprintf('%d-%d', y0, y1);
      return
   end
   w = leg_windows.(model);
   if isstruct(w) && isfield(w, 'start')
      tag = sprintf('%s..%s', ...
         char(string(w.start, 'yyyy-MM-dd')), ...
         char(string(w.end, 'yyyy-MM-dd')));
   else
      tag = char(string(w));   % e.g. "skipped: <reason>"
   end
end

function s = compactYears(years)
   %COMPACTYEARS Render a sorted year list as compact contiguous ranges.
   if isempty(years)
      s = "none";
      return
   end
   years = sort(reshape(years, 1, []));
   parts = strings(0, 1);
   run_start = years(1);
   prev = years(1);
   for k = 2:numel(years)
      if years(k) == prev + 1
         prev = years(k);
         continue
      end
      parts(end+1) = rangeTag(run_start, prev); %#ok<AGROW>
      run_start = years(k);
      prev = years(k);
   end
   parts(end+1) = rangeTag(run_start, prev);
   s = strjoin(parts, ', ');
end

function t = rangeTag(a, b)
   if a == b
      t = string(a);
   else
      t = sprintf('%d-%d', a, b);
   end
end
