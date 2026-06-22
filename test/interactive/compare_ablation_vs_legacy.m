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
   %  (demo/data/eval/kan{l,m}_ablation_daily.mat). The NEW eval builder now
   %  fills the GOF-reference role that bead icemodel-dz2.2 originally assigned
   %  to those staged artifacts, so this check quantifies how well the new
   %  builder reproduces the legacy reference.
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
   %                 (default demo/data/eval under the repo root)
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
      options.legacy_root = fullfile(repo_root, "demo", "data", "eval");
   end
   if options.figdir == ""
      options.figdir = fullfile(repo_root, "test", "interactive", "figures");
   end
   if options.savefigs && ~isfolder(options.figdir)
      mkdir(options.figdir);
   end

   results = struct();
   for s = options.sites
      results.(s) = compareSite(s, options);
   end
   reportResults(results);
end

%% ----------------------------------------------------------------- per site
function r = compareSite(site, options)
   %COMPARESITE NEW buildPromiceData ablation vs the legacy daily artifact.

   r = struct("skipped", true, "reason", "", "stats", struct(), ...
      "figure", "");

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
   la = legacy_a(tf);
   na = new_a(loc(tf));
   gap = new_gap(loc(tf));
   if numel(common_t) < 30
      r.reason = sprintf("%s: only %d overlapping days", site, numel(common_t));
      return
   end

   % --- rebaseline both to the common-window start (relative to install). ---
   la = la - firstFinite(la);
   na = na - firstFinite(na);

   % --- difference stats on direct-observation days only. ---
   valid = ~isnan(la) & ~isnan(na) & ~gap;
   d = na(valid) - la(valid);                % new minus legacy

   r.skipped = false;
   r.stats.window_start = common_t(1);
   r.stats.window_end = common_t(end);
   r.stats.n_days = numel(common_t);
   r.stats.n_valid = nnz(valid);
   r.stats.n_gap_excluded = nnz(gap & ~isnan(la) & ~isnan(na));
   r.stats.bias = mean(d, "omitnan");
   r.stats.rmse = sqrt(mean(d.^2, "omitnan"));
   r.stats.corr = corr(na(valid), la(valid), "rows", "complete");
   r.stats.total_ablation_new = na(find(valid, 1, "last"));
   r.stats.total_ablation_legacy = la(find(valid, 1, "last"));
   r.stats.total_ablation_diff = ...
      r.stats.total_ablation_new - r.stats.total_ablation_legacy;
   [r.stats.max_abs_diff, imax] = max(abs(d));
   vt = common_t(valid);
   r.stats.max_abs_diff_date = vt(imax);

   if options.savefigs
      r.figure = makeFigure(site, common_t, la, na, gap, r.stats, options.figdir);
   end
end

%% ------------------------------------------------------------------ figure
function figpath = makeFigure(site, t, la, na, gap, stats, figdir)
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
   title(sprintf("%s ablation: new vs legacy (rebaselined to window start)", ...
      upper(site)));

   subplot(2, 1, 2);
   d = na - la;
   plot(t, d, "-", "LineWidth", 1.0); grid on; hold on
   yline(stats.bias, "--", "bias", "LabelHorizontalAlignment", "left");
   ylabel("new - legacy [m]"); xlabel("time");
   ttl = sprintf("bias=%.3f m  RMSE=%.3f m  r=%.4f  total_diff=%.3f m  n=%d (gap-excl %d)", ...
      stats.bias, stats.rmse, stats.corr, stats.total_ablation_diff, ...
      stats.n_valid, stats.n_gap_excluded);
   title(ttl, "Interpreter", "none");

   figpath = fullfile(figdir, sprintf("ablation_vs_legacy_%s.png", site));
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

function root = repoRoot()
   %REPOROOT Repo root from this file: test/interactive/<file>.m -> up two.
   root = fileparts(fileparts(fileparts(mfilename("fullpath"))));
end

function reportResults(results)
   %REPORTRESULTS One-line summary per site.
   fprintf("\n--- new PROMICE ablation vs legacy staged ablation ---\n");
   fn = fieldnames(results);
   for k = 1:numel(fn)
      r = results.(fn{k});
      if r.skipped
         fprintf("  %-5s SKIP  (%s)\n", fn{k}, r.reason);
      else
         s = r.stats;
         msg = sprintf("  %-5s n=%d (gap-excl %d)  bias=%.3f m  RMSE=%.3f m  r=%.4f  total_new=%.3f legacy=%.3f diff=%.3f m\n", ...
            fn{k}, s.n_valid, s.n_gap_excluded, s.bias, s.rmse, s.corr, ...
            s.total_ablation_new, s.total_ablation_legacy, ...
            s.total_ablation_diff);
         fprintf("%s", msg);
      end
   end
   fprintf("\n");
end
