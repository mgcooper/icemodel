function results = compare_forcing_vs_legacy(options)
   %COMPARE_FORCING_VS_LEGACY New forcing builders vs legacy ak4 artifacts.
   %
   %  results = compare_forcing_vs_legacy()
   %  results = compare_forcing_vs_legacy(builders=["mar"])
   %  results = compare_forcing_vs_legacy(legacy_root="/path/to/runoff/data")
   %
   %  User-facing comparison of the current icemodel.forcing builders against the
   %  pre-refactor LEGACY ak4 forcing artifacts staged under the runoff project.
   %  This is NOT part of the formal unit suite: it depends on legacy artifacts
   %  (runoff/data/icemodel/...) and S03 raw-source archives that are not durable
   %  fixtures, so it lives in test/interactive and is run by hand.
   %
   %  The checks are STATISTICAL agreement (correlation + bias), not cell-exact
   %  equality, because the legacy cell selection (ncrowcol nearest-match on a
   %  curvilinear grid; catchment-interpolated point for MERRA) cannot be
   %  reproduced cell-for-cell by the current single-cell builders. The thresholds
   %  here mirror the gates the formal forcing tests carried before the
   %  legacy comparison was moved out (owning ExecPlan 2026-06-12 forcing-builder).
   %
   %  Each builder returns a struct with the per-channel correlation / bias and a
   %  logical pass flag; a builder whose legacy artifact or raw source is missing
   %  is reported as skipped rather than failed.
   %
   %  Options:
   %    builders     which comparisons to run (default ["mar","merra"])
   %    legacy_root  root of the legacy artifact tree
   %                 (default /Users/mattcooper/MATLAB/projects/runoff/data)
   %
   % See also: icemodel.forcing.buildMarMet, icemodel.forcing.buildMerraData,
   %  test/unit/test_forcing_mar, test/unit/test_forcing_merra

   arguments
      options.builders (1,:) string = ["mar", "merra"]
      options.legacy_root (1,1) string = ...
         "/Users/mattcooper/MATLAB/projects/runoff/data"
   end

   results = struct();
   for b = options.builders
      switch b
         case "mar"
            results.mar = compareMar(options.legacy_root);
         case "merra"
            results.merra = compareMerra(options.legacy_root);
         otherwise
            warning("compare_forcing_vs_legacy:unknownBuilder", ...
               "unknown builder '%s' (expected mar|merra)", b);
      end
   end

   reportResults(results);
end

%% ---------------------------------------------------------------------- MAR
function r = compareMar(legacy_root)
   %COMPAREMAR buildMarMet vs the legacy ak4 MAR artifact.
   %
   % The legacy ak4 artifacts cannot be reproduced cell-exactly: their cell
   % selection came from ncrowcol's independent row/column nearest-match on a
   % curvilinear grid (stored Lat/Lon metadata are mutually inconsistent), and no
   % single cell of the current MAR archive reproduces the stored series (best
   % whole-year match deviates 6.4 K in tair). The gate is therefore statistical
   % agreement at the one-cell-offset scale.

   r = struct("skipped", true, "passed", false, "reason", "", "stats", struct());

   source_dir = locateMarSource();
   if source_dir == ""
      r.reason = "MAR source data not available (S03 unmounted/spun down, no cache)";
      return
   end
   legacy_file = fullfile(legacy_root, ...
      "icemodel/input/met/met_ak4_mar_2009_1hr.mat");
   if ~isfile(legacy_file)
      r.reason = "legacy ak4 artifact not found: " + legacy_file;
      return
   end

   met = icemodel.forcing.buildMarMet([67.1556, -49.9226], 2009, ...
      source_dir=source_dir);
   legacy = load(legacy_file, "met").("met");
   [tf, loc] = alignTimes(met.Time, legacy.Time);

   r.skipped = false;
   r.stats.n_aligned = sum(tf);
   r.stats.r_tair = corr(met.tair(tf), legacy.tair(loc(tf)), "rows", "complete");
   r.stats.r_lwd  = corr(met.lwd(tf),  legacy.lwd(loc(tf)),  "rows", "complete");
   r.stats.r_psfc = corr(met.psfc(tf), legacy.psfc(loc(tf)), "rows", "complete");
   r.stats.bias_tair = mean(met.tair(tf) - legacy.tair(loc(tf)), "omitnan");
   r.stats.bias_swd  = mean(met.swd(tf)  - legacy.swd(loc(tf)),  "omitnan");

   r.passed = r.stats.n_aligned == 8760 ...
      && r.stats.r_tair > 0.98 ...
      && r.stats.r_lwd > 0.95 ...
      && r.stats.r_psfc > 0.95 ...
      && abs(r.stats.bias_tair) < 2 ...
      && abs(r.stats.bias_swd) < 5;
end

%% -------------------------------------------------------------------- MERRA
function r = compareMerra(legacy_root)
   %COMPAREMERRA buildMerraData vs the legacy ak4_merra artifact.
   %
   % Nearest cell vs the legacy catchment-interpolated point; rh additionally
   % carries the vapor-kernel change; swd carries the SWGDN-vs-SWGNT decision
   % (the legacy artifact derived swd as SWGNT/(1-SNICEALB), inflating it). The
   % gate is correlation per channel plus a small tair bias.

   r = struct("skipped", true, "passed", false, "reason", "", "stats", struct());

   source_dir = locateMerraSource();
   if source_dir == ""
      r.reason = "MERRA-2 source data not available (S03 unmounted/spun down, no cache)";
      return
   end
   legacy_file = fullfile(legacy_root, ...
      "icemodel/input/userdata/ak4_merra_2013.mat");
   if ~isfile(legacy_file)
      r.reason = "legacy ak4 artifact not found: " + legacy_file;
      return
   end

   Data = icemodel.forcing.buildMerraData([67.1556, -49.9226], 2013, ...
      source_dir=source_dir);
   legacy = load(legacy_file, "Data").("Data");
   [tf, loc] = alignTimes(Data.Time, legacy.Time);

   rr = @(v) corr(Data.(v)(tf), legacy.(v)(loc(tf)), "rows", "complete");

   r.skipped = false;
   r.stats.n_aligned = sum(tf);
   r.stats.r_tair   = rr("tair");
   r.stats.r_psfc   = rr("psfc");
   r.stats.r_runoff = rr("runoff");
   r.stats.r_ppt    = rr("ppt");
   r.stats.r_albedo = rr("albedo");
   r.stats.r_swd    = rr("swd");
   r.stats.r_rh     = rr("rh");
   r.stats.bias_tair = mean(Data.tair(tf) - legacy.tair(loc(tf)), "omitnan");

   r.passed = r.stats.n_aligned == 8760 ...
      && r.stats.r_tair > 0.999 ...
      && r.stats.r_psfc > 0.999 ...
      && r.stats.r_runoff > 0.99 ...
      && r.stats.r_ppt > 0.98 ...
      && r.stats.r_albedo > 0.99 ...
      && r.stats.r_swd > 0.95 ...
      && r.stats.r_rh > 0.75 ...
      && abs(r.stats.bias_tair) < 1;
end

%% ------------------------------------------------------------------ helpers
function source_dir = locateMarSource()
   %LOCATEMARSOURCE First readable MAR source dir, or "" when none.
   candidates = [ ...
      "/Volumes/S03/DATA/greenland/mar3p11/RUH2", ...
      string(fullfile(icemodel.getpath("data"), "forcing", "mar"))];
   hasdata = arrayfun(@(p) isfolder(p) ...
      && ~isempty(dir(fullfile(p, "MARv3.11*.nc"))), candidates);
   hit = candidates(hasdata);
   source_dir = "";
   if ~isempty(hit), source_dir = hit(1); end
end

function source_dir = locateMerraSource()
   %LOCATEMERRASOURCE First readable MERRA-2 source dir, or "" when none.
   candidates = [ ...
      "/Volumes/S03/DATA/merra2/1hrly/ncfiles", ...
      string(fullfile(icemodel.getpath("data"), "forcing", "merra2"))];
   hasdata = arrayfun(@(p) ~isempty(dir(fullfile(p, "slv", "*_Nx.*.nc4*"))), ...
      candidates);
   hit = candidates(hasdata);
   source_dir = "";
   if ~isempty(hit), source_dir = hit(1); end
end

function [tf, loc] = alignTimes(t_new, t_legacy)
   %ALIGNTIMES Index-align the new (UTC) series onto the naive legacy axis.
   t_new.TimeZone = "";
   [tf, loc] = ismember(t_new, t_legacy);
end

function reportResults(results)
   %REPORTRESULTS Print a one-line PASS/SKIP/FAIL summary per builder.
   fprintf("\n--- new forcing builders vs legacy ak4 artifacts ---\n");
   fn = fieldnames(results);
   for k = 1:numel(fn)
      r = results.(fn{k});
      if r.skipped
         fprintf("  %-6s SKIP  (%s)\n", fn{k}, r.reason);
      elseif r.passed
         fprintf("  %-6s PASS  (n=%d aligned)\n", fn{k}, r.stats.n_aligned);
      else
         fprintf("  %-6s FAIL  (n=%d aligned; inspect results.%s.stats)\n", ...
            fn{k}, r.stats.n_aligned, fn{k});
      end
   end
   fprintf("\n");
end
