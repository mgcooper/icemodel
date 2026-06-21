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
   % over the newly-migrated pypromice Level-3 product (the QC'd L3 surface
   % and ice-temperature channels, read not derived):
   %
   %   (a) Snow depth [m] from the L3 snow_height channel. Should be >= 0,
   %       seasonal (winter peak, late-summer ~0), and O(0.1-1 m). Persistent
   %       negatives or multi-metre spikes are suspect.
   %
   %   (b) Cumulative ablation [m, positive down] from the L3 ice-surface
   %       height (z_ice_surf), or z_surf_combined for accumulation-zone
   %       stations that ship no z_ice_surf. The source channel is annotated
   %       on the panel and recorded in the sanity table. Ablation should
   %       step up each melt season and stay flat through winter.
   %
   %   (c) The subsurface ice-temperature string tice1..ticeN [degC] as a
   %       line set (shallow -> deep). Shallow sensors swing with the seasonal
   %       cycle; deeper sensors damp and lag; all sit below ~0 degC after the
   %       dictionary [-80, 1] C clamp applied in readPromiceAws.
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

      % --- (a) Snow depth (L3 snow_height) --------------------------------
      ax1 = nexttile(tl);
      if hasChan(Data, "snow_depth")
         plot(ax1, t, Data.snow_depth, '-', 'LineWidth', 1.0)
      end
      yline(ax1, 0, ':', 'Color', [0.5 0.5 0.5])
      ylabel(ax1, 'snow depth [m]')
      title(ax1, '(a) snow depth (L3 snow\_height)')
      grid(ax1, 'on')

      % --- (b) Cumulative ablation (L3 surface height) --------------------
      ax2 = nexttile(tl);
      if hasChan(Data, "ablation")
         plot(ax2, t, Data.ablation, '-', 'LineWidth', 1.0)
      end
      yline(ax2, 0, ':', 'Color', [0.5 0.5 0.5])
      ylabel(ax2, 'cumulative ablation [m, +down]')
      title(ax2, sprintf('(b) cumulative ablation  [source: %s]', ...
         meta.ablation_source), 'Interpreter', 'none')
      grid(ax2, 'on')

      % --- (c) Subsurface ice-temperature string --------------------------
      ax3 = nexttile(tl);
      have = tice_names(ismember(tice_names, ...
         string(Data.Properties.VariableNames)));
      hold(ax3, 'on')
      for k = 1:numel(have)
         % Kelvin -> degC for readability (channels stored in kelvin).
         plot(ax3, t, Data.(have{k}) - 273.15, '-', 'LineWidth', 0.7)
      end
      yline(ax3, 0, ':', 'Color', [0.5 0.5 0.5])
      hold(ax3, 'off')
      ylabel(ax3, 'ice temperature [degC]')
      xlabel(ax3, 'time (UTC)')
      title(ax3, '(c) subsurface temperature string tice1..ticeN (shallow->deep)')
      if ~isempty(have)
         legend(ax3, have, 'Location', 'eastoutside', 'Interpreter', 'none')
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
   sd = colFinite(Data, "snow_depth");
   ab = colFinite(Data, "ablation");
   have = tice_names(ismember(tice_names, ...
      string(Data.Properties.VariableNames)));

   % Snow depth metrics.
   if isempty(sd)
      sd_med = NaN; sd_max = NaN; sd_neg = NaN;
   else
      sd_med = median(sd); sd_max = max(sd);
      sd_neg = nnz(sd < -0.02);   % tolerate sub-cm noise around zero
   end

   % Ablation metrics: total surface lowering and a coarse monotonicity
   % fraction (share of steps that do not move upward, i.e. surface lowering
   % or flat). Accumulation-dominated sites legitimately fail monotonicity.
   if isempty(ab)
      ab_total = NaN; ab_mono = NaN;
   else
      ab_total = ab(end) - ab(1);
      d = diff(ab);
      ab_mono = nnz(d >= -1e-6) / max(numel(d), 1);
   end

   % tice metrics after the [-80, 1] C clamp.
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
      string(upper(site)), ...
      string(Data.Time(1)), string(Data.Time(end)), ...
      round(span_days), height(Data), ...
      round(sd_med, 2), round(sd_max, 2), sd_neg, ...
      round(ab_total, 2), round(ab_mono, 2), string(meta.ablation_source), ...
      numel(have), round(ti_min, 1), round(ti_max, 1), round(100 * ti_warm, 1), ...
      "", ...
      'VariableNames', {'station', 'start', 'stop', 'span_d', 'nrows', ...
      'sd_med', 'sd_max', 'sd_neg', 'abl_total', 'abl_mono', 'abl_source', ...
      'n_tice', 'tice_min', 'tice_max', 'tice_warmpct', 'flags'});

   row.flags = buildFlags(row, frequency);
end

function row = failedRow(site, ME)
   %FAILEDROW Placeholder row for a station the builder could not build.
   nanr = @() NaN;
   row = table( ...
      string(upper(site)), missing, missing, nanr(), nanr(), ...
      nanr(), nanr(), nanr(), nanr(), nanr(), "(build failed)", ...
      nanr(), nanr(), nanr(), nanr(), ...
      string("FAIL: " + regexprep(ME.message, '\s+', ' ')), ...
      'VariableNames', {'station', 'start', 'stop', 'span_d', 'nrows', ...
      'sd_med', 'sd_max', 'sd_neg', 'abl_total', 'abl_mono', 'abl_source', ...
      'n_tice', 'tice_min', 'tice_max', 'tice_warmpct', 'flags'});
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

   % Snow depth. Sub-cm negatives are tolerated in sd_neg; only flag when a
   % nontrivial number of samples sit persistently below the surface.
   if isnan(row.sd_med)
      f(end+1) = "NO_SNOWDEPTH";
   elseif row.sd_neg > 50
      f(end+1) = "SNOW_NEG";
   end
   if ~isnan(row.sd_max) && row.sd_max > 5
      f(end+1) = "SNOW_HIGH";
   end

   % Ablation.
   if isnan(row.abl_total)
      f(end+1) = "NO_ABLATION";
   end
   if row.abl_source == "L3 z_surf_combined"
      f(end+1) = "ABL_ZSURF";          % accumulation site, fallback source
   elseif row.abl_source == "none (no L3 surface-height channel)"
      f(end+1) = "ABL_MISSING";
   end
   % Net accumulation (negative total ablation) is fine on accumulation
   % sites but worth flagging for the user's eye.
   if ~isnan(row.abl_total) && row.abl_total < -0.5
      f(end+1) = "ABL_NET_ACCUM";
   end

   % tice string. The clamp ceiling is +1 C, so isolated near-melt touches at
   % the ceiling are normal; flag only when warm samples are a meaningful
   % SHARE of the string (a persistently-warm/biased thermistor record).
   if row.n_tice == 0
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
   fprintf(['  %-8s span %5sd  n=%-6s  sd[med %5s max %5s neg %4s]  ' ...
      'abl[tot %6s src %-18s]  tice[%2s %5s..%5sC warm%%%5s]  %s\n'], ...
      row.station, num(row.span_d), num(row.nrows), ...
      num(row.sd_med), num(row.sd_max), num(row.sd_neg), ...
      num(row.abl_total), char(row.abl_source), num(row.n_tice), ...
      num(row.tice_min), num(row.tice_max), num(row.tice_warmpct), row.flags);
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
      "NO_SNOWDEPTH",     "no finite snow_depth (L3 snow_height) samples"
      "SNOW_NEG",         "> 50 snow_depth samples persistently < -2 cm"
      "SNOW_HIGH",        "snow_depth max > 5 m (suspect)"
      "NO_ABLATION",      "no finite ablation samples"
      "ABL_ZSURF",        "ablation from z_surf_combined (accumulation site)"
      "ABL_MISSING",      "no L3 surface-height channel for ablation"
      "ABL_NET_ACCUM",    "net negative ablation (surface raised > 0.5 m)"
      "NO_TICE",          "no subsurface temperature sensors"
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
