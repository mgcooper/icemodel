function plot_thf_cases(options)
   %PLOT_THF_CASES Load saved THF validation results and produce figures.
   %
   %  plot_thf_cases()
   %  plot_thf_cases(group_idx=1:2)
   %  plot_thf_cases(datafile='/path/to/file.mat')
   %  plot_thf_cases(datafile='/path/to/file.mat', group_idx=3)
   %
   % Loads results saved by run_thf_cases and produces two figures per
   % site/year group (four groups total):
   %
   %   Figure 1 — 2×2 scatter (Tsfc, NetR, SHF, LHF):
   %     x = hourly observations, y = icemodel; both schemes overlaid with
   %     a 1:1 reference line, a best-fit line, and R² in the legend.
   %
   %   Figure 2 — icemodel.plot.enbal timeseries:
   %     Daily 7-day moving-mean; obs = solid, scheme A = dashed,
   %     scheme B = dotted; one color per variable (van As palette).
   %
   % Both figures are saved as PNG to test/interactive/figures/.
   %
   % group_idx selects which site/year groups to plot (default: all).
   % Groups are sorted alphabetically by 'site_simyears' key, e.g.:
   %   1: kanl_2015-2016
   %   2: kanl_2016
   %   3: kanm_2015-2016
   %   4: kanm_2016

   arguments
      options.datafile (1,:) char = ''
      options.group_idx (1,:) double = []
      options.save_figs (1,:) logical = true
   end

   % Default data file resolved at runtime (mfilename not valid in arguments block).
   if isempty(options.datafile)
      options.datafile = fullfile(fileparts(mfilename('fullpath')), 'data', ...
         'thf_validation_results.mat');
   end

   loaded = load(options.datafile, 'cases', 'run_data');
   cases = loaded.cases;
   run_data = loaded.run_data;

   figdir = fullfile(fileparts(mfilename('fullpath')), 'figures');
   if ~exist(figdir, 'dir')
      mkdir(figdir);
   end

   % Group cases by site + simulation year span (e.g. 'kanm_2016').
   keys = arrayfun(@(c) ...
      sprintf('%s_%s', c.site, strjoin(string(c.simyears), '-')), cases, ...
      'UniformOutput', false);
   groups = unique(keys);

   % Apply group_idx filter if provided.
   if ~isempty(options.group_idx)
      groups = groups(options.group_idx);
      fprintf('Plotting groups: %s\n', strjoin(groups, ', '));
   end

   for g = 1:numel(groups)
      gdata = run_data(strcmp(keys, groups{g}));
      gname = groups{g};
      gtitle = strrep(gname, '_', ' ');

      % --- Figure 1: scatter comparison ---
      fig_sc = plot_scatter_comparison(gdata, gtitle);

      if options.save_figs
         saveas(fig_sc, fullfile(figdir, [gname '_scatter.png']));
      end

      % --- Figure 2: enbal timeseries comparison ---
      % Pass grouptitle so enbal renders a 2-line title without sgtitle overlap.
      [ice1a, ice1b, met_obs] = extract_scheme_pair(gdata);
      icemodel.plot.enbal(ice1a, met_obs, ice1b, ...
         'labels', {'Bulk-Richardson', 'Monin-Obukhov'}, ...
         'grouptitle', gtitle);

      if options.save_figs
         saveas(gcf, fullfile(figdir, [gname '_timeseries.png']));
      end
   end
end

% =========================================================================
% Scatter comparison
% =========================================================================

function fig = plot_scatter_comparison(run_data, gtitle)
   %PLOT_SCATTER_COMPARISON 2×2 scatter panel: both schemes vs observations.
   %
   % x = hourly observations, y = icemodel. For each variable and scheme: a
   % point cloud, a least-squares best-fit line, and R² in the legend.
   % A 1:1 dashed reference line spans the full data range.

   import icemodel.helpers.rmttleapinds

   vars = {'tsfc', 'netr', 'shf', 'lhf'};
   panel_titles = {'T_{sfc} (°C)', 'Net Radiation (W m^{-2})', ...
      'SHF (W m^{-2})', 'LHF (W m^{-2})'};

   % Per-scheme colors (teal / amber — avoids blue and red).
   cmap.bulk_richardson = struct('color', [0.07 0.58 0.55], ...
      'name', 'Bulk-Richardson');
   cmap.monin_obukhov = struct('color', [0.85 0.50 0.10], ...
      'name', 'Monin-Obukhov');

   n_runs = numel(run_data);

   fig = figure('Name', [gtitle ' — scatter'], 'Position', [100 100 1400 900]);
   tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
   title(tl, gtitle, 'FontWeight', 'bold');

   for vi = 1:numel(vars)
      vn = vars{vi};
      ax = nexttile(tl);
      hold(ax, 'on');

      obs_series = cell(n_runs, 1);
      model_series = cell(n_runs, 1);
      colors = zeros(n_runs, 3);
      names = strings(n_runs, 1);
      n_valid = 0;

      for irun = 1:n_runs
         d = run_data{irun};
         pr = cmap.(d.c.scheme);

         TT = sync_hourly(d.ice1, d.met, vn);
         if isempty(TT)
            continue
         end
         keep = isfinite(TT.model) & isfinite(TT.ref);
         obs = TT.ref(keep);
         mdl = TT.model(keep);
         if numel(obs) < 2
            continue
         end
         n_valid = n_valid + 1;
         obs_series{n_valid} = obs;
         model_series{n_valid} = mdl;
         colors(n_valid, :) = pr.color;
         names(n_valid) = pr.name;
      end

      if n_valid > 0
         icemodel.plot.scatterplot( ...
            obs_series(1:n_valid), model_series(1:n_valid), ...
            axes=ax, ...
            display_name=names(1:n_valid), ...
            colors=colors(1:n_valid, :), ...
            x_label="Observations", ...
            y_label="icemodel", ...
            marker_size=18, ...
            marker_face_alpha=0.30, ...
            line_width=2.0, ...
            legend_location="northwest");
      end

      title(ax, panel_titles{vi});
   end
end

% =========================================================================
% Helpers
% =========================================================================

function [ice1a, ice1b, met_obs] = extract_scheme_pair(gdata)
   %EXTRACT_SCHEME_PAIR Pull out the two scheme runs and their shared met.

   ice1a = [];
   ice1b = [];
   met_obs = [];
   for di = 1:numel(gdata)
      d = gdata{di};
      if isempty(met_obs)
         met_obs = d.met;
      end
      if strcmp(d.c.scheme, 'bulk_richardson')
         ice1a = d.ice1;
      else
         ice1b = d.ice1;
      end
   end
end

function TT = sync_hourly(ice1, met, varname)
   %SYNC_HOURLY Synchronize one model and one met variable at hourly resolution.
   %
   % Returns empty if the variable is absent from either timetable.

   import icemodel.helpers.rmttleapinds

   if ~ismember(varname, ice1.Properties.VariableNames) || ...
         ~ismember(varname, met.Properties.VariableNames)
      TT = [];
      return
   end

   % Retime to hourly (data is typically already hourly; rmttleapinds removes
   % Feb 29 rows so leap years do not produce mismatched time bases).
   mdl_h = rmttleapinds(retime(ice1, 'hourly', 'mean'));
   ref_h = rmttleapinds(retime(met, 'hourly', 'mean'));

   mdl_tt = timetable(mdl_h.Time, mdl_h.(varname), 'VariableNames', {'model'});
   ref_tt = timetable(ref_h.Time, ref_h.(varname), 'VariableNames', {'ref'});
   TT = synchronize(mdl_tt, ref_tt, 'intersection');
end
