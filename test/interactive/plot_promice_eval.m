function plot_promice_eval(options)
   %PLOT_PROMICE_EVAL Sanity-check the PROMICE evaluation channels.
   %
   %  plot_promice_eval()
   %  plot_promice_eval(sites=["kanm" "kanl"])
   %  plot_promice_eval(source_dir="/Volumes/S03/DATA/greenland/geus/aws/v3")
   %  plot_promice_eval(save_figs=false)
   %
   % Builds the PROMICE evaluation Data timetable
   % (icemodel.forcing.buildPromiceData) for each requested site and
   % produces one diagnostic figure per site with three stacked panels:
   %
   %   (a) Snow depth [m] derived from the sonic boom height, overlaid on
   %       the raw boom height so the September-median reference / clamp
   %       derivation in buildPromiceData/snowDepthFromBoom is visually
   %       checkable. Snow depth should be >= 0, seasonal (winter peak,
   %       late-summer ~0), and O(0.1-1 m).
   %
   %   (b) Cumulative ablation [m, positive down] from the despiked
   %       pressure-transducer record, overlaid on the raw transducer
   %       depth so the despike / reset-removal in cumulativeAblation is
   %       visually checkable. Ablation should be monotonic-ish upward
   %       (surface lowering accumulates), step up each melt season, and
   %       stay roughly flat through winter.
   %
   %   (c) The subsurface ice-temperature string tice1..tice8 [degC] as a
   %       line set (shallow -> deep). Shallow strings should swing with
   %       the seasonal cycle; deeper strings damp and lag; all should sit
   %       below ~0 degC firn/ice temperatures.
   %
   % This is a verification/diagnostics tool for the firn evaluation-data
   % work: the PROMICE eval channels (snow depth, ablation, the tice
   % string) were never independently sanity-checked before they became
   % load-bearing for the PROMICE expansion. It is not a unit test; the
   % automated contract lives in test/unit/test_forcing_promice.m.
   %
   % Figures are written as PNG to test/interactive/figures/ (gitignored).
   %
   % Name-value
   %  sites      : station ids to plot (default ["kanm" "kanl"])
   %  source_dir : PROMICE v3 NetCDF dir; "" lets the builder resolve it
   %               (S03 external drive layout or the local cache)
   %  frequency  : "daily" (default, readable) or "hourly"
   %  save_figs  : write PNGs to test/interactive/figures/ (default true)
   %
   % See also: icemodel.forcing.buildPromiceData,
   %  icemodel.forcing.readPromiceAws, test_forcing_promice

   arguments
      options.sites (1, :) string = ["kanm", "kanl"]
      options.source_dir (1, 1) string = ""
      options.frequency (1, 1) string ...
         {mustBeMember(options.frequency, ["hourly", "daily"])} = "daily"
      options.save_figs (1, 1) logical = true
   end

   % Resolve the figure output directory relative to this file so the
   % script is location-independent (mfilename not valid in arguments).
   figdir = fullfile(fileparts(mfilename('fullpath')), 'figures');
   if options.save_figs && ~isfolder(figdir)
      mkdir(figdir)
   end

   % The tice string, shallow to deep, used by panel (c).
   tice_names = arrayfun(@(n) sprintf('tice%d', n), 1:8, ...
      'UniformOutput', false);

   for site = options.sites
      % Build the evaluation Data timetable for this site. Skip the site
      % (rather than abort the whole sweep) if its source file is absent.
      try
         [Data, meta] = icemodel.forcing.buildPromiceData(site, ...
            source_dir=options.source_dir, frequency=options.frequency);
      catch ME
         fprintf('SKIP %-6s: %s\n', site, ...
            regexprep(ME.message, '\s+', ' '));
         continue
      end

      % Print a compact numeric sanity summary alongside the figure so the
      % eyeball check is backed by magnitudes/signs in the console log.
      printSummary(site, Data, meta, tice_names);

      t = Data.Time;
      fig = figure('Name', sprintf('PROMICE eval - %s', upper(site)), ...
         'Position', [100, 100, 900, 900]);
      tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', ...
         'Padding', 'compact');
      title(tl, sprintf('PROMICE evaluation channels - %s (%s)', ...
         upper(site), options.frequency), 'Interpreter', 'none');

      % --- (a) Snow depth vs raw boom height -------------------------------
      ax1 = nexttile(tl);
      yyaxis(ax1, 'left')
      plot(ax1, t, Data.snow_depth, '-', 'LineWidth', 1.2)
      ylabel(ax1, 'snow depth [m]')
      yyaxis(ax1, 'right')
      if ismember("boom_height", string(Data.Properties.VariableNames))
         plot(ax1, t, Data.boom_height, '-', 'LineWidth', 0.6)
      end
      ylabel(ax1, 'boom height [m]')
      title(ax1, '(a) snow depth (derived) vs raw boom height')
      grid(ax1, 'on')

      % --- (b) Cumulative ablation vs raw transducer depth -----------------
      ax2 = nexttile(tl);
      yyaxis(ax2, 'left')
      plot(ax2, t, Data.ablation, '-', 'LineWidth', 1.2)
      ylabel(ax2, 'cumulative ablation [m, +down]')
      yyaxis(ax2, 'right')
      if ismember("transducer_depth", string(Data.Properties.VariableNames))
         plot(ax2, t, Data.transducer_depth, '-', 'LineWidth', 0.6)
      end
      ylabel(ax2, 'transducer depth [m]')
      title(ax2, '(b) cumulative ablation (despiked) vs raw transducer depth')
      grid(ax2, 'on')

      % --- (c) Subsurface ice-temperature string ---------------------------
      ax3 = nexttile(tl);
      have = tice_names(ismember(tice_names, ...
         string(Data.Properties.VariableNames)));
      hold(ax3, 'on')
      for k = 1:numel(have)
         % Convert kelvin -> degC for readability; the channels are stored
         % in kelvin per the builder docstring.
         plot(ax3, t, Data.(have{k}) - 273.15, '-', 'LineWidth', 0.8)
      end
      hold(ax3, 'off')
      ylabel(ax3, 'ice temperature [degC]')
      xlabel(ax3, 'time (UTC)')
      title(ax3, '(c) subsurface temperature string tice1..tice8 (shallow->deep)')
      if ~isempty(have)
         legend(ax3, have, 'Location', 'eastoutside', 'Interpreter', 'none')
      end
      grid(ax3, 'on')

      % Save a small PNG for the committed verification record.
      if options.save_figs
         figfile = fullfile(figdir, ...
            sprintf('promice_eval_%s.png', lower(erase(site, "_"))));
         exportgraphics(fig, figfile, 'Resolution', 150)
         fprintf('wrote %s\n', figfile);
      end
   end
end

%% Local functions
function printSummary(site, Data, meta, tice_names)
   %PRINTSUMMARY Console sanity summary backing the visual eyeball check.
   sd = Data.snow_depth(isfinite(Data.snow_depth));
   ab = Data.ablation(isfinite(Data.ablation));

   fprintf('\n=== %s  (%s, %d rows) ===\n', upper(site), ...
      meta.station_name, height(Data));
   fprintf('  window      : %s -> %s\n', ...
      string(Data.Time(1)), string(Data.Time(end)));

   % Snow depth: expect non-negative, sub-meter-to-few-meter seasonal.
   if isempty(sd)
      fprintf('  snow_depth  : (all missing)\n');
   else
      fprintf('  snow_depth  : min %.2f  median %.2f  max %.2f m  (n=%d)\n', ...
         min(sd), median(sd), max(sd), numel(sd));
   end

   % Ablation: expect positive total surface lowering of O(meters/season).
   if isempty(ab)
      fprintf('  ablation    : (all missing)\n');
   else
      fprintf('  ablation    : min %.2f  max %.2f m  (n=%d)\n', ...
         min(ab), max(ab), numel(ab));
   end

   % tice string: report degC range per available depth.
   have = tice_names(ismember(tice_names, ...
      string(Data.Properties.VariableNames)));
   for k = 1:numel(have)
      v = Data.(have{k})(isfinite(Data.(have{k}))) - 273.15;
      if ~isempty(v)
         fprintf('  %-7s     : min %6.1f  median %6.1f  max %6.1f degC\n', ...
            have{k}, min(v), median(v), max(v));
      end
   end
end
