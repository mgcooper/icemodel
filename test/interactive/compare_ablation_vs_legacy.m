function results = compare_ablation_vs_legacy(options)
   %COMPARE_ABLATION_VS_LEGACY New PROMICE ablation vs legacy staged ablation.
   %
   %  results = compare_ablation_vs_legacy()
   %  results = compare_ablation_vs_legacy(sites=["kanl"])
   %  results = compare_ablation_vs_legacy(savefigs=false)
   %
   %  One-off, USER-FACING comparison of the NEW PROMICE ablation built by
   %  icemodel.forcing.buildPromiceData (reads the L3 z_ice_surf channel with
   %  site-type branching + gap flags) against the LEGACY pre-refactor staged
   %  ablation artifacts at KAN_L and KAN_M
   %  (test/assets/legacy_ablation/kan{l,m}_ablation_daily.mat). The NEW eval
   %  builder now fills the GOF-reference role that bead icemodel-dz2.2
   %  originally assigned to those staged artifacts, so this check quantifies
   %  how well the new builder reproduces the legacy reference, and whether a
   %  historical model-data evaluation would change if the obs reference were
   %  swapped from the legacy series to the new one.
   %
   %  Two comparisons are produced for each site:
   %    1. FULL-RECORD overlap (the legacy 2009-2018 span), as before.
   %    2. Targeted MELT-SEASON WINDOWS (2015/2016, Jun-Oct and Jul-Oct), the
   %       windows over which historical KAN_L/KAN_M model evaluations were run.
   %       Each window rebaselines BOTH series to their value at the window
   %       start t1 (the legacy plotPromice(...,'refstart',t1) convention:
   %       subtract the t1 value so both series start at 0 at t1), then reports
   %       how much a GOF metric would shift if the obs reference changed.
   %
   %  This is NOT part of the formal unit suite. It depends on the committed
   %  PROMICE L3 cache (data/verification/promice/day) and the staged legacy
   %  artifacts, and writes figures to the gitignored test/interactive/figures/.
   %
   %  Both series are CUMULATIVE surface lowering relative to installation, so
   %  each is rebaselined to its own value at the common-window start before
   %  comparison (the absolute offset is meaningless; only the shape and the
   %  total ablation over the window are physical).
   %
   %  Options:
   %    sites     which sites to compare (default ["kanl","kanm"])
   %    legacy_root  dir holding the legacy *_ablation_daily.mat artifacts
   %                 (default test/assets/legacy_ablation under the repo root)
   %    figdir    output figure dir (default test/interactive/figures)
   %    savefigs  write png overlays (default true)
   %
   % See also: icemodel.forcing.buildPromiceData,
   %  test/interactive/compare_forcing_vs_legacy

   arguments
      options.sites (1, :) string = ["kanl", "kanm"]
      options.legacy_root (1, 1) string = ""
      options.figdir (1, 1) string = ""
      options.savefigs (1, 1) logical = true
   end

   repo_root = repoRoot();
   if options.legacy_root == ""
      options.legacy_root = fullfile(repo_root, "test", "assets", ...
         "legacy_ablation");
   end
   if options.figdir == ""
      options.figdir = fullfile(repo_root, "test", "interactive", "figures");
   end
   if options.savefigs && ~isfolder(options.figdir)
      mkdir(options.figdir);
   end

   % Targeted melt-season evaluation windows (naive, year/month/day). Each is
   % rebaselined to the first valid sample at or after the window start.
   windows = meltSeasonWindows();

   results = struct();
   for s = options.sites
      results.(s) = compareSite(s, windows, options);
   end
   reportResults(results);
end

%% ----------------------------------------------------------------- windows
function w = meltSeasonWindows()
   %MELTSEASONWINDOWS Targeted historical-evaluation windows (naive datetimes).
   w = struct( ...
      "name", {"2015_jun_oct", "2015_jul_oct", "2016_jun_oct", "2016_jul_oct"}, ...
      "t1",   {datetime(2015,6,1), datetime(2015,7,1), ...
               datetime(2016,6,1), datetime(2016,7,1)}, ...
      "t2",   {datetime(2015,10,1), datetime(2015,10,1), ...
               datetime(2016,10,1), datetime(2016,10,1)});
end

%% ----------------------------------------------------------------- per site
function r = compareSite(site, windows, options)
   %COMPARESITE NEW buildPromiceData ablation vs the legacy daily artifact.

   r = struct("skipped", true, "reason", "", "stats", struct(), ...
      "figure", "", "windows", struct([]));

   % --- legacy staged ablation (cumulative, daily, positive down). ---
   legacy_file = fullfile(options.legacy_root, site + "_ablation_daily.mat");
   if ~isfile(legacy_file)
      r.reason = "legacy artifact not found: " + legacy_file;
      return
   end
   L = load(legacy_file);
   legacy = L.ablation;                      % timetable, var 'ablation'
   legacy_t = legacy.Properties.RowTimes;
   legacy_t.TimeZone = "";                    % naive axis for daily alignment
   legacy_a = legacy.ablation;

   % --- NEW PROMICE ablation from the L3 z_ice_surf channel. ---
   try
      Data = icemodel.forcing.buildPromiceData(siteId(site), frequency="daily");
   catch err
      r.reason = "buildPromiceData(" + site + ") failed: " + err.message;
      return
   end
   if ~ismember("ablation", string(Data.Properties.VariableNames))
      r.reason = site + " is not an ablation site (no z_ice_surf / ablation)";
      return
   end
   new_t = Data.Properties.RowTimes;
   new_t.TimeZone = "";                      % drop UTC tz for naive alignment
   new_a = Data.ablation;

   % Gap-bridged samples are not direct observations; exclude them from the
   % stats so divergence reflects real disagreement, not interpolation.
   new_gap = false(size(new_a));
   if ismember("surface_height_flag", string(Data.Properties.VariableNames))
      new_gap = Data.surface_height_flag > 0;
   end

   % --- common daily axis (align on calendar date; daily means may be
   % timestamped at a different hour than the legacy midnight stamp). ---
   legacy_day = dateshift(legacy_t, "start", "day");
   new_day = dateshift(new_t, "start", "day");
   [tf, loc] = ismember(legacy_day, new_day);
   common_t = legacy_t(tf);
   common_day = legacy_day(tf);
   la = legacy_a(tf);
   na = new_a(loc(tf));
   gap = new_gap(loc(tf));
   if numel(common_t) < 30
      r.reason = sprintf("%s: only %d overlapping days", site, numel(common_t));
      return
   end

   r.skipped = false;

   % --- (1) FULL-RECORD comparison (rebaselined to the overlap start). ---
   r.stats = windowStats(common_t, ...
      la - firstFinite(la), na - firstFinite(na), gap);

   if options.savefigs
      r.figure = makeFigure(site, common_t, ...
         la - firstFinite(la), na - firstFinite(na), gap, r.stats, ...
         options.figdir, "full", "full overlap (2009-2018)");
   end

   % --- (2) targeted melt-season windows (rebaselined to window start t1). --
   for k = 1:numel(windows)
      sel = common_day >= windows(k).t1 & common_day <= windows(k).t2;
      wt = common_t(sel);
      lw = la(sel);
      nw = na(sel);
      gw = gap(sel);
      ws = windowStats(wt, lw - firstFinite(lw), nw - firstFinite(nw), gw);
      ws.name = windows(k).name;
      ws.t1 = windows(k).t1;
      ws.t2 = windows(k).t2;
      if isempty(r.windows)
         r.windows = ws;
      else
         r.windows(end+1) = ws;
      end
      if options.savefigs && ws.n_valid >= 5
         makeFigure(site, wt, lw - firstFinite(lw), nw - firstFinite(nw), ...
            gw, ws, options.figdir, windows(k).name, ...
            sprintf("%s -> %s", string(windows(k).t1, "dd-MMM-uuuu"), ...
            string(windows(k).t2, "dd-MMM-uuuu")));
      end
   end
end

%% ----------------------------------------------------------- window stats
function s = windowStats(t, la, na, gap)
   %WINDOWSTATS Agreement stats for one rebaselined window (new vs legacy).
   %  la, na are already rebaselined (start at 0 at window start). The stats
   %  use direct-observation days only (drop NaN and gap-bridged samples).
   s = struct();
   valid = ~isnan(la) & ~isnan(na) & ~gap;
   d = na(valid) - la(valid);                % new minus legacy

   s.window_start = firstOr(t, NaT);
   s.window_end = lastOr(t, NaT);
   s.n_days = numel(t);
   s.n_valid = nnz(valid);
   s.n_gap_excluded = nnz(gap & ~isnan(la) & ~isnan(na));
   if s.n_valid == 0
      [s.bias, s.rmse, s.corr] = deal(NaN);
      [s.total_ablation_new, s.total_ablation_legacy] = deal(NaN);
      [s.total_ablation_diff, s.total_ablation_relpct] = deal(NaN);
      [s.max_abs_diff, s.max_abs_diff_date] = deal(NaN, NaT);
      return
   end
   s.bias = mean(d, "omitnan");
   s.rmse = sqrt(mean(d.^2, "omitnan"));
   if s.n_valid >= 2
      s.corr = corr(na(valid), la(valid), "rows", "complete");
   else
      s.corr = NaN;
   end
   s.total_ablation_new = na(find(valid, 1, "last"));
   s.total_ablation_legacy = la(find(valid, 1, "last"));
   s.total_ablation_diff = s.total_ablation_new - s.total_ablation_legacy;
   % Relative shift of the window total ablation reference (new vs legacy):
   % this is how much a normalized GOF (e.g. NSE/relative bias) would move if
   % the obs reference were swapped from legacy to new over this window.
   if abs(s.total_ablation_legacy) > eps
      s.total_ablation_relpct = ...
         100 * s.total_ablation_diff / s.total_ablation_legacy;
   else
      s.total_ablation_relpct = NaN;
   end
   [s.max_abs_diff, imax] = max(abs(d));
   vt = t(valid);
   s.max_abs_diff_date = vt(imax);
end

%% ------------------------------------------------------------------ figure
function figpath = makeFigure(site, t, la, na, gap, stats, figdir, tag, span_lbl)
   %MAKEFIGURE Overlay (top) + difference (bottom); save png.
   f = figure("Visible", "off", "Position", [100 100 1000 700]);

   subplot(2, 1, 1);
   plot(t, la, "-", "LineWidth", 1.4, "DisplayName", "legacy staged"); hold on
   plot(t, na, "-", "LineWidth", 1.0, "DisplayName", "new buildPromiceData");
   if any(gap)
      plot(t(gap), na(gap), ".", "MarkerSize", 4, ...
         "DisplayName", "gap-bridged (excluded)");
   end
   grid on; legend("Location", "northwest");
   ylabel("cumulative ablation [m, +down]");
   title(sprintf("%s ablation: new vs legacy, rebaselined to start (%s)", ...
      upper(site), span_lbl), "Interpreter", "none");

   subplot(2, 1, 2);
   d = na - la;
   plot(t, d, "-", "LineWidth", 1.0); grid on; hold on
   yline(stats.bias, "--", "bias", "LabelHorizontalAlignment", "left");
   ylabel("new - legacy [m]"); xlabel("time");
   ttl = sprintf(['bias=%.3f m  RMSE=%.3f m  r=%.4f  total new=%.3f ' ...
      'legacy=%.3f diff=%.3f m (%.1f%%)  n=%d (gap-excl %d)'], ...
      stats.bias, stats.rmse, stats.corr, stats.total_ablation_new, ...
      stats.total_ablation_legacy, stats.total_ablation_diff, ...
      relpct(stats), stats.n_valid, stats.n_gap_excluded);
   title(ttl, "Interpreter", "none");

   figpath = fullfile(figdir, sprintf("ablation_vs_legacy_%s_%s.png", ...
      site, tag));
   exportgraphics(f, figpath, "Resolution", 150);
   close(f);
end

%% ------------------------------------------------------------------ helpers
function id = siteId(site)
   %SITEID Map a compact alias (kanl) to the PROMICE station id (KAN_L).
   id = upper(site);
   id = replace(id, "KANL", "KAN_L");
   id = replace(id, "KANM", "KAN_M");
   id = replace(id, "KANU", "KAN_U");
end

function v0 = firstFinite(v)
   %FIRSTFINITE First finite value of a vector (0 when none).
   idx = find(isfinite(v), 1);
   if isempty(idx)
      v0 = 0;
   else
      v0 = v(idx);
   end
end

function v = firstOr(x, default)
   %FIRSTOR First element of x, or default when x is empty.
   if isempty(x); v = default; else; v = x(1); end
end

function v = lastOr(x, default)
   %LASTOR Last element of x, or default when x is empty.
   if isempty(x); v = default; else; v = x(end); end
end

function p = relpct(stats)
   %RELPCT Relative total-ablation difference percent (NaN-safe field access).
   if isfield(stats, "total_ablation_relpct")
      p = stats.total_ablation_relpct;
   else
      p = NaN;
   end
end

function root = repoRoot()
   %REPOROOT Repo root from this file: test/interactive/<file>.m -> up two.
   root = fileparts(fileparts(fileparts(mfilename("fullpath"))));
end

function reportResults(results)
   %REPORTRESULTS Full-record summary + per-window table per site.
   fprintf("\n--- new PROMICE ablation vs legacy staged ablation ---\n");
   fn = fieldnames(results);
   for k = 1:numel(fn)
      r = results.(fn{k});
      if r.skipped
         fprintf("  %-5s SKIP  (%s)\n", fn{k}, r.reason);
      else
         s = r.stats;
         fprintf(['  %-5s FULL  n=%d (gap-excl %d)  bias=%.3f m  RMSE=%.3f m' ...
            '  r=%.4f  total new=%.3f legacy=%.3f diff=%.3f m (%.1f%%)\n'], ...
            fn{k}, s.n_valid, s.n_gap_excluded, s.bias, s.rmse, s.corr, ...
            s.total_ablation_new, s.total_ablation_legacy, ...
            s.total_ablation_diff, relpct(s));
      end
   end

   % Per-window table (the historical-evaluation question).
   fprintf("\n--- melt-season windows (rebaselined to window start) ---\n");
   fprintf("%-6s %-13s %5s %5s %7s %7s %7s %8s %8s %8s %7s\n", ...
      "site", "window", "n", "gapx", "bias", "RMSE", "r", ...
      "tot_new", "tot_leg", "diff", "rel%");
   for k = 1:numel(fn)
      r = results.(fn{k});
      if r.skipped || isempty(r.windows); continue; end
      for j = 1:numel(r.windows)
         w = r.windows(j);
         fprintf(['%-6s %-13s %5d %5d %7.3f %7.3f %7.4f %8.3f %8.3f ' ...
            '%8.3f %7.1f\n'], fn{k}, w.name, w.n_valid, w.n_gap_excluded, ...
            w.bias, w.rmse, w.corr, w.total_ablation_new, ...
            w.total_ablation_legacy, w.total_ablation_diff, ...
            relpct(w));
      end
   end
   fprintf("\n");
end
