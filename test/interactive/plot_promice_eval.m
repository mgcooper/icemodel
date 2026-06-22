function summary = plot_promice_eval(options)
   %PLOT_PROMICE_EVAL Sanity-check the QC'd PROMICE evaluation channels.
   %
   %  plot_promice_eval()                       % all staged stations
   %  plot_promice_eval(sites="all")            % explicit all
   %  plot_promice_eval(sites=["KAN_M" "KAN_L"])
   %  plot_promice_eval(save_figs=false)
   %  T = plot_promice_eval();                  % return the summary table
   %
   % Builds the PROMICE evaluation Data timetable
   % (icemodel.forcing.buildPromiceData) for each requested station and
   % produces one diagnostic figure per station with three stacked panels
   % over the pypromice Level-3 product (the QC'd L3 surface and
   % ice-temperature channels, read not derived). The surface panels are
   % SITE-TYPE aware (buildPromiceData branches on z_ice_surf presence):
   %
   %   (a) ABLATION sites: snow depth [m] from the clamped L3 snow_height
   %       channel (>= 0, seasonal, O(0.1-1 m)). ACCUMULATION sites: this
   %       panel is the net surface_height [m, positive up] from
   %       z_surf_combined, since they ship no true snow-over-ice depth.
   %
   %   (b) ABLATION sites: cumulative ablation [m, positive down] from
   %       z_ice_surf (step up each melt season, flat through winter).
   %       ACCUMULATION sites: the same surface_height series with the
   %       GAP-BRIDGED segments (surface_height_flag==1) overplotted in red so
   %       interpolated, non-observational stretches are visually distinct and
   %       can be excluded from model comparison. The source channel is
   %       annotated and recorded in the sanity table.
   %
   %   (c) Subsurface temperature: the standardized 10 m channel tice10m
   %       [degC] as a heavy line (the PRIMARY evaluation channel) overlaid on
   %       the depth-tagged thermistor string tice1..ticeN [degC] (shallow ->
   %       deep), all below ~0 degC after the dictionary [-80, 1] C clamp and
   %       the surfaced-thermistor discard in readPromiceAws.
   %
   % This is a verification/diagnostics tool for the firn evaluation-data
   % work: the migrated L3 eval channels (snow depth, ablation, the tice
   % string) need a per-station eyeball check before they become load-bearing
   % for model development. It is not a unit test; the automated contract
   % lives in test/unit/test_forcing_promice.m.
   %
   % Alongside the figures, a per-station SANITY SUMMARY TABLE is printed to
   % the console and saved as a markdown file in the (gitignored) figures
   % directory, flagging stations whose channels look wrong or sparse so the
   % user can prioritise.
   %
   % Figures (PNG) and the summary markdown are written to
   % test/interactive/figures/ (gitignored).
   %
   % Name-value
   %  sites      : station ids to plot, or "all" (default) to sweep every
   %               hour/<STATION>_hour.nc under source_dir. Accepts canonical
   %               ids ("KAN_M") or compact aliases ("kanm").
   %  source_dir : PROMICE L3 NetCDF dir or its parent holding hour/. Defaults
   %               to the staged product data/verification/promice (NO
   %               external-drive fallback).
   %  frequency  : "daily" (default, readable) or "hourly".
   %  save_figs  : write PNGs + summary to figures/ (default true).
   %
   % Output
   %  summary    : table, one row per station, with the sanity metrics (also
   %               printed and, when save_figs, saved as markdown).
   %
   % See also: icemodel.forcing.buildPromiceData,
   %  icemodel.forcing.readPromiceAws, test_forcing_promice

   arguments
      options.sites (1, :) string = "all"
      options.source_dir (1, 1) string = ""
      options.frequency (1, 1) string ...
         {mustBeMember(options.frequency, ["hourly", "daily"])} = "daily"
      options.save_figs (1, 1) logical = true
   end

   % Default source_dir to the staged verification product (no /Volumes
   % fallback): resolve it relative to the repo so the builders read the
   % migrated L3 bundle directly.
   source_dir = options.source_dir;
   if source_dir == ""
      source_dir = string(fullfile(icemodel.internal.fullpath('data'), ...
         'verification', 'promice'));
   end

   % Expand sites="all" to every station file staged under hour/.
   sites = options.sites;
   if isscalar(sites) && lower(sites) == "all"
      sites = discoverStations(source_dir);
      if isempty(sites)
         error('plot_promice_eval:noStations', ...
            'no hour/*.nc station files found under %s', source_dir)
      end
   end

   % Resolve the figure output directory relative to this file so the script
   % is location-independent (mfilename not valid in arguments).
   figdir = fullfile(fileparts(mfilename('fullpath')), 'figures');
   if options.save_figs && ~isfolder(figdir)
      mkdir(figdir)
   end

   % The tice string, shallow to deep, used by panel (c). The L3 string runs
   % to ~11 sensors; allow headroom.
   tice_names = arrayfun(@(n) sprintf('tice%d', n), 1:12, ...
      'UniformOutput', false);

   rows = cell(numel(sites), 1);
   for s = 1:numel(sites)
      site = sites(s);

      % Build the evaluation Data timetable for this station. Record the
      % failure (rather than abort the whole sweep) if the file is absent or
      % the builder errors, so the summary reflects every station attempted.
      try
         [Data, meta] = icemodel.forcing.buildPromiceData(site, ...
            source_dir=source_dir, frequency=options.frequency);
      catch ME
         fprintf('FAIL %-8s: %s\n', site, regexprep(ME.message, '\s+', ' '));
         rows{s} = failedRow(site, ME);
         continue
      end

      % Collect the per-station sanity metrics (drives the summary table and
      % the console print below).
      rows{s} = sanityRow(site, Data, meta, tice_names, options.frequency);
      printSummary(rows{s});

      t = Data.Time;
      fig = figure('Name', sprintf('PROMICE eval - %s', upper(site)), ...
         'Position', [100, 100, 900, 950], 'Visible', 'off');
      tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', ...
         'Padding', 'compact');
      title(tl, sprintf('PROMICE evaluation channels - %s (%s)', ...
         upper(site), options.frequency), 'Interpreter', 'none');

      is_ablation = isfield(meta, 'site_surface_type') ...
         && string(meta.site_surface_type) == "ablation";

      % --- (a) snow depth (ablation) or surface height (accumulation) -----
      ax1 = nexttile(tl);
      if is_ablation && hasChan(Data, "snow_depth")
         plot(ax1, t, Data.snow_depth, '-', 'LineWidth', 1.0)
         ylabel(ax1, 'snow depth [m]')
         title(ax1, '(a) snow depth (ablation site, L3 snow\_height, clamped >= 0)')
      elseif hasChan(Data, "surface_height")
         plot(ax1, t, Data.surface_height, '-', 'LineWidth', 1.0)
         ylabel(ax1, 'surface height [m, +up]')
         title(ax1, '(a) net surface height (accumulation site, L3 z\_surf\_combined)')
      end
      yline(ax1, 0, ':', 'Color', [0.5 0.5 0.5])
      grid(ax1, 'on')

      % --- (b) surface-height channel with gap-bridged segments flagged ---
      ax2 = nexttile(tl);
      if is_ablation && hasChan(Data, "ablation")
         surf = Data.ablation; surflbl = 'cumulative ablation [m, +down]';
      elseif hasChan(Data, "surface_height")
         surf = Data.surface_height; surflbl = 'surface height [m, +up]';
      else
         surf = nan(height(Data), 1); surflbl = 'surface height [m]';
      end
      hold(ax2, 'on')
      plot(ax2, t, surf, '-', 'Color', [0 0.3 0.7], 'LineWidth', 1.0)
      % Overplot gap-bridged samples (flag==1) in red so non-observational
      % segments are visually distinct and excludable from model comparison.
      if hasChan(Data, "surface_height_flag")
         gap = Data.surface_height_flag == 1;
         sg = surf; sg(~gap) = NaN;
         plot(ax2, t, sg, '.', 'Color', [0.85 0.1 0.1], 'MarkerSize', 4)
         legend(ax2, {'observed', 'gap-bridged (flag=1)'}, ...
            'Location', 'best', 'Interpreter', 'none')
      end
      yline(ax2, 0, ':', 'Color', [0.5 0.5 0.5])
      hold(ax2, 'off')
      ylabel(ax2, surflbl)
      title(ax2, sprintf('(b) surface height  [source: %s]', ...
         meta.surface_source), 'Interpreter', 'none')
      grid(ax2, 'on')

      % --- (c) tice10m (primary) + depth-tagged thermistor string ---------
      ax3 = nexttile(tl);
      have = tice_names(ismember(tice_names, ...
         string(Data.Properties.VariableNames)));
      hold(ax3, 'on')
      for k = 1:numel(have)
         % Kelvin -> degC for readability (channels stored in kelvin).
         plot(ax3, t, Data.(have{k}) - 273.15, '-', 'LineWidth', 0.6, ...
            'Color', [0.6 0.6 0.6])
      end
      legentries = have;
      if hasChan(Data, "tice10m")
         plot(ax3, t, Data.tice10m - 273.15, '-', 'Color', [0 0 0], ...
            'LineWidth', 1.6)
         legentries = [{'tice10m (PRIMARY)'}, have];
      end
      yline(ax3, 0, ':', 'Color', [0.5 0.5 0.5])
      hold(ax3, 'off')
      ylabel(ax3, 'ice temperature [degC]')
      xlabel(ax3, 'time (UTC)')
      title(ax3, ['(c) subsurface temperature: tice10m (primary) + ' ...
         'depth-tagged string (shallow->deep)'])
      if ~isempty(legentries)
         legend(ax3, legentries, 'Location', 'eastoutside', ...
            'Interpreter', 'none')
      end
      grid(ax3, 'on')

      if options.save_figs
         figfile = fullfile(figdir, ...
            sprintf('promice_eval_%s.png', lower(erase(site, "_"))));
         exportgraphics(fig, figfile, 'Resolution', 150)
      end
      close(fig)
   end

   summary = vertcat(rows{:});
   summary = sortrows(summary, 'station');

   % Print the table and save the markdown sanity summary.
   printTable(summary);
   if options.save_figs
      writeMarkdown(summary, fullfile(figdir, 'promice_eval_summary.md'), ...
         options.frequency, source_dir);
      fprintf('\nwrote figures + summary to %s\n', figdir);
   end
end

%% Local functions ------------------------------------------------------------
function sites = discoverStations(source_dir)
   %DISCOVERSTATIONS Station ids from every hour/<STATION>_hour.nc file.
   hourdir = fullfile(source_dir, 'hour');
   if ~isfolder(hourdir)
      hourdir = source_dir;
   end
   files = dir(fullfile(hourdir, '*_hour.nc'));
   names = string({files.name});
   sites = erase(names, "_hour.nc");
   sites = sort(sites);
end

function tf = hasChan(Data, name)
   %HASCHAN True if the timetable carries (any finite sample of) a channel.
   tf = ismember(name, string(Data.Properties.VariableNames)) ...
      && any(isfinite(Data.(char(name))));
end

function row = sanityRow(site, Data, meta, tice_names, frequency)
   %SANITYROW One-row sanity table for a built station.
   site_type = "unknown";
   if isfield(meta, 'site_surface_type')
      site_type = string(meta.site_surface_type);
   end
   surf_source = "n/a";
   if isfield(meta, 'surface_source')
      surf_source = string(meta.surface_source);
   end

   % Snow depth metrics (ablation sites only; accumulation sites carry none).
   sd = colFinite(Data, "snow_depth");
   if isempty(sd)
      sd_med = NaN; sd_max = NaN; sd_neg = NaN;
   else
      sd_med = median(sd); sd_max = max(sd);
      sd_neg = nnz(sd < -0.02);   % tolerate sub-cm noise around zero
   end

   % Surface metrics on the site-type channel: ablation (positive down) for
   % ablation sites, surface_height (positive up) for accumulation sites. The
   % "mono" fraction is the share of steps not moving against the expected
   % sign (lowering for ablation; rising for accumulation).
   if site_type == "ablation"
      surfv = colFinite(Data, "ablation");
   else
      surfv = colFinite(Data, "surface_height");
   end
   if isempty(surfv)
      ab_total = NaN; ab_mono = NaN;
   else
      ab_total = surfv(end) - surfv(1);
      d = diff(surfv);
      ab_mono = nnz(d >= -1e-6) / max(numel(d), 1);
   end

   % Gap-bridged sample share on the surface channel.
   gf = colFinite(Data, "surface_height_flag");
   if isempty(gf)
      gap_pct = NaN;
   else
      gap_pct = 100 * nnz(gf == 1) / numel(gf);
   end

   % tice metrics after the [-80, 1] C clamp and surfaced-thermistor discard.
   have = tice_names(ismember(tice_names, ...
      string(Data.Properties.VariableNames)));
   has_tice10m = ismember("tice10m", string(Data.Properties.VariableNames)) ...
      && any(isfinite(Data.tice10m));
   tv = [];
   for k = 1:numel(have)
      tv = [tv; colFinite(Data, have{k}) - 273.15]; %#ok<AGROW>
   end
   if isempty(tv)
      ti_min = NaN; ti_max = NaN; ti_warm = NaN;
   else
      ti_min = min(tv); ti_max = max(tv);
      % Warm-sample share: fraction of finite tice samples above 0.5 C. The
      % dictionary clamp ceiling is +1 C, so brief near-melt touches read at
      % the ceiling and are normal; a large warm SHARE is the suspect signal.
      ti_warm = nnz(tv > 0.5) / numel(tv);
   end

   span_days = days(Data.Time(end) - Data.Time(1));

   row = table( ...
      string(upper(site)), site_type, ...
      string(Data.Time(1)), string(Data.Time(end)), ...
      round(span_days), height(Data), ...
      round(sd_med, 2), round(sd_max, 2), sd_neg, ...
      round(ab_total, 2), round(ab_mono, 2), surf_source, round(gap_pct, 1), ...
      numel(have), has_tice10m, round(ti_min, 1), round(ti_max, 1), ...
      round(100 * ti_warm, 1), "", ...
      'VariableNames', {'station', 'site_type', 'start', 'stop', 'span_d', ...
      'nrows', 'sd_med', 'sd_max', 'sd_neg', 'surf_total', 'surf_mono', ...
      'surf_source', 'gap_pct', 'n_tice', 'has_tice10m', 'tice_min', ...
      'tice_max', 'tice_warmpct', 'flags'});

   row.flags = buildFlags(row, frequency);
end

function row = failedRow(site, ME)
   %FAILEDROW Placeholder row for a station the builder could not build.
   nanr = @() NaN;
   row = table( ...
      string(upper(site)), "unknown", missing, missing, nanr(), nanr(), ...
      nanr(), nanr(), nanr(), nanr(), nanr(), "(build failed)", nanr(), ...
      nanr(), false, nanr(), nanr(), nanr(), ...
      string("FAIL: " + regexprep(ME.message, '\s+', ' ')), ...
      'VariableNames', {'station', 'site_type', 'start', 'stop', 'span_d', ...
      'nrows', 'sd_med', 'sd_max', 'sd_neg', 'surf_total', 'surf_mono', ...
      'surf_source', 'gap_pct', 'n_tice', 'has_tice10m', 'tice_min', ...
      'tice_max', 'tice_warmpct', 'flags'});
end

function flags = buildFlags(row, frequency)
   %BUILDFLAGS Concatenate the suspect-channel flags for one station.
   f = string.empty;

   % Short / sparse record (< 1 yr, or sparse coverage vs nominal cadence).
   if row.span_d < 365
      f(end+1) = "SHORT(<1yr)";
   end
   step_h = 24; if frequency == "hourly"; step_h = 1; end
   coverage = row.nrows / max((row.span_d * 24 / step_h), 1);
   if coverage < 0.5
      f(end+1) = "SPARSE";
   end

   % Snow depth (ablation sites only; NO_SNOWDEPTH on accumulation sites is
   % EXPECTED, not a defect, so it is not flagged there).
   if row.site_type == "ablation"
      if isnan(row.sd_med)
         f(end+1) = "NO_SNOWDEPTH";
      elseif row.sd_neg > 50
         f(end+1) = "SNOW_NEG";
      end
      if ~isnan(row.sd_max) && row.sd_max > 5
         f(end+1) = "SNOW_HIGH";
      end
   end

   % Surface-height channel.
   if isnan(row.surf_total)
      f(end+1) = "NO_SURFACE";
   end
   if row.surf_source == "none (no L3 surface-height channel)"
      f(end+1) = "SURF_MISSING";
   end
   % A large gap-bridged share means much of the surface series is not a
   % direct observation (exclude those segments before model comparison).
   if ~isnan(row.gap_pct) && row.gap_pct > 25
      f(end+1) = "GAP_HEAVY(>25%)";
   end
   % Net ablation sign sanity: an ablation site that NET rises, or an
   % accumulation site that NET lowers, is worth the user's eye.
   if row.site_type == "ablation" && ~isnan(row.surf_total) ...
         && row.surf_total < -0.5
      f(end+1) = "ABL_NET_RISE";
   elseif row.site_type == "accumulation" && ~isnan(row.surf_total) ...
         && row.surf_total < -0.5
      f(end+1) = "ACC_NET_LOWER";
   end

   % tice string. The clamp ceiling is +1 C, so isolated near-melt touches at
   % the ceiling are normal; flag only when warm samples are a meaningful
   % SHARE of the string (a persistently-warm/biased thermistor record).
   if row.n_tice == 0 && ~row.has_tice10m
      f(end+1) = "NO_TICE";
   end
   if ~isnan(row.tice_warmpct) && row.tice_warmpct > 2
      f(end+1) = "TICE_WARM(>2%>0.5C)";
   end

   if isempty(f)
      flags = "ok";
   else
      flags = strjoin(f, ",");
   end
end

function v = colFinite(Data, name)
   %COLFINITE Finite values of a channel (empty if the channel is absent).
   name = char(name);
   if ismember(name, Data.Properties.VariableNames)
      v = Data.(name);
      v = v(isfinite(v));
   else
      v = [];
   end
end

function printSummary(row)
   %PRINTSUMMARY Compact per-station console line backing the eyeball check.
   fprintf(['  %-8s %-12s span %5sd  n=%-6s  sd[med %5s neg %4s]  ' ...
      'surf[tot %6s gap%%%5s src %-18s]  tice[%2s %5s..%5sC warm%%%5s]  %s\n'], ...
      row.station, char(row.site_type), num(row.span_d), num(row.nrows), ...
      num(row.sd_med), num(row.sd_neg), ...
      num(row.surf_total), num(row.gap_pct), char(row.surf_source), ...
      num(row.n_tice), num(row.tice_min), num(row.tice_max), ...
      num(row.tice_warmpct), row.flags);
end

function printTable(summary)
   %PRINTTABLE Echo the full sanity table + the flagged subset to console.
   fprintf('\n===== PROMICE eval sanity summary (%d stations) =====\n', ...
      height(summary));
   disp(summary);
   flagged = summary(summary.flags ~= "ok", :);
   fprintf('----- FLAGGED stations (%d) -----\n', height(flagged));
   for i = 1:height(flagged)
      fprintf('  %-8s : %s\n', flagged.station(i), flagged.flags(i));
   end
end

function writeMarkdown(summary, mdfile, frequency, source_dir)
   %WRITEMARKDOWN Save the sanity table as a small markdown file.
   fid = fopen(mdfile, 'w');
   if fid < 0
      warning('plot_promice_eval:write', 'cannot write %s', mdfile);
      return
   end
   cleanup = onCleanup(@() fclose(fid));
   fprintf(fid, '# PROMICE L3 evaluation-channel sanity summary\n\n');
   fprintf(fid, '- generated: %s\n', string(datetime("now")));
   fprintf(fid, '- frequency: %s\n', frequency);
   fprintf(fid, '- source_dir: %s\n', source_dir);
   fprintf(fid, '- stations: %d\n\n', height(summary));

   cols = summary.Properties.VariableNames;
   fprintf(fid, '| %s |\n', strjoin(cols, ' | '));
   fprintf(fid, '|%s|\n', repmat(' --- |', 1, numel(cols)));
   for i = 1:height(summary)
      cells = strings(1, numel(cols));
      for c = 1:numel(cols)
         cells(c) = cellString(summary.(cols{c})(i));
      end
      fprintf(fid, '| %s |\n', strjoin(cells, ' | '));
   end

   flagged = summary(summary.flags ~= "ok", :);
   fprintf(fid, '\n## Flagged stations (%d)\n\n', height(flagged));
   for i = 1:height(flagged)
      fprintf(fid, '- **%s**: %s\n', flagged.station(i), flagged.flags(i));
   end

   fprintf(fid, '\n## Flag legend\n\n');
   legend = [
      "SHORT(<1yr)",      "record spans < 1 year"
      "SPARSE",           "row count < 50% of nominal cadence coverage"
      "NO_SNOWDEPTH",     "ablation site with no finite snow_depth samples"
      "SNOW_NEG",         "> 50 snow_depth samples persistently < -2 cm"
      "SNOW_HIGH",        "snow_depth max > 5 m (suspect)"
      "NO_SURFACE",       "no finite surface-height samples (ablation or surface_height)"
      "SURF_MISSING",     "no L3 surface-height channel at all"
      "GAP_HEAVY(>25%)",  "> 25% of surface-height samples gap-bridged (flag=1)"
      "ABL_NET_RISE",     "ablation site net-rose > 0.5 m (station transition / data issue)"
      "ACC_NET_LOWER",    "accumulation site net-lowered > 0.5 m (suspect)"
      "NO_TICE",          "no subsurface temperature (no thermistor string and no tice10m)"
      "TICE_WARM(>2%>0.5C)", "> 2% of tice samples above 0.5 C (warm bias)"];
   for i = 1:size(legend, 1)
      fprintf(fid, '- `%s`: %s\n', legend(i, 1), legend(i, 2));
   end
end

function s = cellString(v)
   %CELLSTRING Markdown-cell string for a scalar table value.
   if isstring(v) || ischar(v)
      s = string(v);
      if ismissing(s)
         s = "-";
      end
   elseif isdatetime(v)
      if ismissing(v)
         s = "-";
      else
         s = string(v, 'yyyy-MM-dd');
      end
   elseif isnumeric(v)
      if isnan(v)
         s = "-";
      else
         s = string(v);
      end
   else
      s = string(v);
   end
end

function s = num(v)
   %NUM Compact string for a possibly-NaN numeric (for fprintf %s fields).
   if isnan(v)
      s = "-";
   else
      s = string(v);
   end
end
