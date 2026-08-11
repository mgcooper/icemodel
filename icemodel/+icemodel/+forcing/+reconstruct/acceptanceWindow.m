function [window, proxy_files] = acceptanceWindow(site, kwargs)
   %ACCEPTANCEWINDOW Per-site forcing-ready policy window from staged proxies.
   %
   %  window = icemodel.forcing.reconstruct.acceptanceWindow( ...
   %     "kanm", location=location)
   %
   % Role
   %  The forcing-ready acceptance-window policy. A reconstruction verdict must
   %  not be penalized by missing proxy TIME coverage, which is a staging
   %  limitation, only by fundamental limitations. The window within which
   %  readiness is judged is therefore the continuous union of the per-source
   %  staged proxy files loaded for the site: the left edge is the earliest
   %  selected sample, and the right edge is the latest. The filename date
   %  tokens must contain those endpoints, but they never widen a partial
   %  boundary day beyond the timetable's exact support. This function
   %  validates every candidate file and rejects a disjoint staged union. For
   %  each source the widest file anchors the selection, and validated siblings
   %  that extend beyond it are also pinned, so reconstruction loads every
   %  staged span and the window reports it.
   %
   %  The window derives from validated met artifacts and their window-stamped
   %  FILENAMES at call time. Restaging therefore redefines the policy with no
   %  refill and no stale stored columns: extending MERRA-2 widens the window
   %  automatically. A new source (for example HIRHAM) joins by a deliberate
   %  code change at the marked setopts extension point (supportedProxySources,
   %  POLICY A12) and by staging its met. The stored readiness ledger stays
   %  ABSOLUTE: it records whether the file is forcing-ready. Consumers apply
   %  this window as a read-time view.
   %
   % Name-value
   %  met_dir : selected native source directory or flat met root. Proxy
   %     discovery uses the shared per-source-first, flat-fallback layout
   %     (default data/input/met/promice under the repo).
   %  location : requested target point used to verify every proxy artifact.
   %  opts : reconstruction options (proxy catalog; default setopts()).
   %
   % Returns
   %  window : 1x2 datetime [start, end] of the continuous staged-proxy span;
   %     [NaT, NaT] when no proxy met is staged for the site.
   %  proxy_files : string column of the selected per-source filenames
   %     contributing to the window, for callers that must pin the read-time
   %     view.
   %
   % See also: icemodel.forcing.reconstruct.setopts,
   %  icemodel.forcing.reconstruct.fillPromiceStation

   arguments
      site (1, 1) string {icemodel.forcing.reconstruct.mustBeStationToken}
      kwargs.met_dir (1, 1) string = ""
      kwargs.location (1, 1) struct
      kwargs.opts (1, 1) struct = icemodel.forcing.reconstruct.setopts()
   end

   met_dir = kwargs.met_dir;
   if met_dir == ""
      met_dir = string(fullfile(icemodel.internal.fullpath, ...
         'data', 'input', 'met', 'promice'));
   end
   [~, met_root] = ...
      icemodel.forcing.reconstruct.selectedDataRoot(met_dir);

   % Union span across every staged window file of every catalog source.
   % Validate each artifact before its filename tokens widen the policy window.
   % A stale or mislabeled narrow file is as unsafe as the widest one.
    catalog = kwargs.opts.proxy_catalog;
    window = NaT(1, 2, 'TimeZone', 'UTC');
   % Per-source selections collect in cells because one source can pin
   % several files: the widest file plus any span extenders (POLICY A6 -
   % staged met must widen the product, never shrink it).
    source_files = cell(numel(catalog), 1);
    source_sample_coverage = cell(numel(catalog), 1);
   % The staged-inventory coverage collects per source in the same way, so
   % the union concatenates once after the catalog loop.
   staged_sample_coverage = cell(numel(catalog), 1);
    for k = 1:numel(catalog)
       hits = [];
        for d = icemodel.forcing.helpers.sourceSearchDirs( ...
              met_root, catalog(k).storage)
          candidates = dir(fullfile(d{1}, sprintf( ...
             'met_%s_%s_*_15m.mat', site, catalog(k).storage)));
          valid_name = false(numel(candidates), 1);
          for h = 1:numel(candidates)
             valid_name(h) = ~isempty(regexp(candidates(h).name, ...
                '_(\d{8})_(\d{8})_15m\.mat$', 'once'));
          end
          if any(valid_name)
             hits = candidates(valid_name);
             break
          end
        end
        valid_hit = false(numel(hits), 1);
        hit_sample_windows = NaT(numel(hits), 2, 'TimeZone', 'UTC');
        for h = 1:numel(hits)
         tokens = regexp(hits(h).name, ...
            '_(\d{8})_(\d{8})_15m\.mat$', 'tokens', 'once');
         if isempty(tokens)
            continue
         end
         [proxy, proxy_file] = ...
            icemodel.forcing.reconstruct.loadWidestTimetable(hits(h));
         if isempty(proxy) || ...
               ~icemodel.forcing.reconstruct.proxyArtifactIdentity( ...
               proxy.Properties.UserData, site, kwargs.location, ...
               catalog(k).storage)
            error(['icemodel:reconstruct:acceptanceWindow:' ...
               'proxyIdentityMismatch'], ...
               ['staged proxy metadata does not identify target %s and ' ...
               'source %s: %s'], site, catalog(k).label, proxy_file);
         end
         t0 = datetime(tokens{1}, 'InputFormat', 'yyyyMMdd', ...
            'TimeZone', 'UTC');
          t1 = datetime(tokens{2}, 'InputFormat', 'yyyyMMdd', ...
             'TimeZone', 'UTC');
          proxy_times = proxy.Properties.RowTimes;
          valid_support = false;
          if numel(proxy_times) >= 2
             endpoint_dates = dateshift(proxy_times([1, end]), ...
                'start', 'day');
             valid_support = isequal(endpoint_dates(:).', [t0, t1]) ...
                && all(diff(proxy_times) == minutes(15)) ...
                && all(mod(minute(proxy_times), 15) == 0 ...
                & second(proxy_times) == 0);
          end
          if ~valid_support
             error(['icemodel:reconstruct:acceptanceWindow:' ...
                'proxyWindowMismatch'], ...
                ['staged proxy timetable endpoint dates, 15-minute ' ...
                'cadence, or UTC quarter-hour grid do not match the ' ...
                'filename window for %s: %s'], ...
                site, proxy_file);
         end
         valid_hit(h) = true;
         hit_sample_windows(h, :) = proxy_times([1, end]);
        end
      % Every validated file of this source counts toward the staged
      % inventory, whether or not the selection below picks it.
      staged_sample_coverage{k} = hit_sample_windows(valid_hit, :);
        if any(valid_hit)
           % The widest file anchors the source. Validated siblings that
           % extend coverage beyond it join the selection, so a staged
           % year is never dropped just because a file with a wider
           % duration exists (POLICY A6: staging more proxy met widens
           % the product). Interior overlap still belongs to the anchor.
           spans = hit_sample_windows(:, 2) - hit_sample_windows(:, 1);
           spans(~valid_hit) = -Inf;
           [~, widest] = max(spans);
           extend = valid_hit ...
              & (hit_sample_windows(:, 1) ...
              < hit_sample_windows(widest, 1) ...
              | hit_sample_windows(:, 2) ...
              > hit_sample_windows(widest, 2));
           extend(widest) = false;
           picked = [widest; find(extend)];
           names = strings(numel(picked), 1);
           for f = 1:numel(picked)
              names(f) = string(fullfile(hits(picked(f)).folder, ...
                 hits(picked(f)).name));
           end
           source_files{k} = names;
           source_sample_coverage{k} = hit_sample_windows(picked, :);
        end
     end
     proxy_files = unique(vertcat(strings(0, 1), ...
        source_files{:}), 'stable');
     selected_sample_coverage = vertcat( ...
        NaT(0, 2, 'TimeZone', 'UTC'), source_sample_coverage{:});
   sample_coverage = vertcat( ...
      NaT(0, 2, 'TimeZone', 'UTC'), staged_sample_coverage{:});
    if isempty(sample_coverage)
       return
    end

     % Reject sub-day holes in the staged inventory and in the exact files that
     % reconstruction selects. Filename dates alone cannot prove continuity.
     continuousWindow(sample_coverage, site, minutes(15));
     window = continuousWindow(selected_sample_coverage, site, minutes(15));
end

function window = continuousWindow(coverage, site, adjacency)
   %CONTINUOUSWINDOW Collapse intervals separated by no more than adjacency.
   coverage = sortrows(coverage, 1);
   covered_through = coverage(1, 2);
   for k = 2:size(coverage, 1)
      if coverage(k, 1) > covered_through + adjacency
         error(['icemodel:reconstruct:acceptanceWindow:' ...
            'proxyCoverageGap'], ...
            'staged proxy windows do not continuously cover target %s', site);
      end
      covered_through = max(covered_through, coverage(k, 2));
   end
   window = [coverage(1, 1), covered_through];
end
