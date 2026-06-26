function leg = resolveLegWindows(models, coverage, window_start, window_end)
   %RESOLVELEGWINDOWS Decouple each gridded RCM leg's window from the met window.
   %
   %  leg = icemodel.verification.setup.resolveLegWindows(models, coverage, ...
   %     window_start, window_end)
   %
   %  Resolves, for each requested gridded source, the calendar window/years to
   %  stage by intersecting the requested met window with the source's on-disk
   %  coverage (probed cheaply by promiceSourceCoverage). This is the FAIL-EARLY
   %  gate: a source with no overlap is marked staged=false WITH A REASON before
   %  any NetCDF is opened, so an empty model/window is skipped without entering
   %  an hours-long build (RR3 feedback #1).
   %
   %  Per-source policy (all three intersect the requested window via capLeg):
   %    * MAR / MERRA : met sources. Window = requested window intersected with
   %      on-disk years. When the requested window is unbounded (NaT - the
   %      all-available default), the source's FULL on-disk coverage is used.
   %    * RACMO       : eval/reference Data only (no met), but still intersected
   %      with the window (8fc): a station whose record lies entirely outside
   %      RACMO's coverage gets a SKIPPED RACMO leg rather than a zero-overlap
   %      file; an unbounded window falls back to RACMO's full coverage.
   %
   %  Inputs
   %    models       : string vector subset of ["mar","merra","racmo"] (other
   %                   entries, e.g. "promice", are ignored - PROMICE is not a
   %                   gridded RCM leg).
   %    coverage     : struct from promiceSourceCoverage (per-source year ranges).
   %    window_start : datetime/NaT requested met window start (NaT = unbounded).
   %    window_end   : datetime/NaT requested met window end (NaT = unbounded).
   %
   %  Returns
   %    leg : struct with a field per requested gridded source; each leg is
   %          struct('staged', logical, 'years', vector, 'start', datetime, ...
   %          'end', datetime, 'reason', string). Consume staged/years/start/end.
   %
   % See also: icemodel.verification.setup.promiceSourceCoverage,
   %  icemodel.verification.setup.stageRcmForcing,
   %  icemodel.verification.setup.importPromiceSites

   leg = struct();

   if ismember("mar", models)
      leg.mar = capLeg(coverage.mar, window_start, window_end, "MAR");
   end
   if ismember("merra", models)
      leg.merra = capLeg(coverage.merra, window_start, window_end, "MERRA-2");
   end
   if ismember("racmo", models)
      % RACMO is intersected with the requested window like MAR/MERRA (8fc):
      % staging RACMO's full 2012-2018 span for a station whose observations are
      % entirely outside it (e.g. a 2022+ record) produces a zero-overlap,
      % unusable Data file. capLeg SKIPS on no overlap and CLIPS on partial. An
      % unbounded (NaT) window still falls back to RACMO's full coverage.
      leg.racmo = capLeg(coverage.racmo, window_start, window_end, "RACMO");
   end
end

%% Local helpers
function L = capLeg(cov, window_start, window_end, label)
   %CAPLEG Met leg: requested window intersected with on-disk years.
   if isempty(cov.years)
      L = struct('staged', false, 'years', [], 'start', NaT, 'end', NaT, ...
         'reason', sprintf('%s absent (%s)', label, cov.reason));
      return
   end

   % An unbounded (NaT) requested window means "all available" - fall back to
   % the source's full on-disk coverage (mirrors ownLeg). This is the SUMup
   % all-available default; PROMICE always passes a bounded station window.
   if isnat(window_start) || isnat(window_end)
      L = ownLeg(cov, label);
      return
   end

   req_years = year(window_start):year(window_end);
   years = intersect(req_years, cov.years);
   if isempty(years)
      L = struct('staged', false, 'years', [], 'start', NaT, 'end', NaT, ...
         'reason', sprintf('%s on-disk %d-%d has no overlap with requested %d-%d', ...
         label, cov.year_min, cov.year_max, req_years(1), req_years(end)));
      return
   end
   y0 = max(year(window_start), min(years));
   y1 = min(year(window_end), max(years));
   t1 = max(window_start, icemodel.verification.setup.ensureUtc( ...
      sprintf('%d-01-01', y0)));
   t2 = min(window_end, icemodel.verification.setup.ensureUtc( ...
      sprintf('%d-12-31 23:00:00', y1)));
   L = struct('staged', true, 'years', years(years >= y0 & years <= y1), ...
      'start', t1, 'end', t2, 'reason', "");
end

function L = ownLeg(cov, label)
   %OWNLEG Full on-disk coverage, decoupled from any requested window.
   if isempty(cov.years)
      L = struct('staged', false, 'years', [], 'start', NaT, 'end', NaT, ...
         'reason', sprintf('%s absent (%s)', label, cov.reason));
      return
   end
   t1 = icemodel.verification.setup.ensureUtc( ...
      sprintf('%d-01-01', cov.year_min));
   t2 = icemodel.verification.setup.ensureUtc( ...
      sprintf('%d-12-31 23:00:00', cov.year_max));
   L = struct('staged', true, 'years', cov.years, 'start', t1, 'end', t2, ...
      'reason', "");
end
