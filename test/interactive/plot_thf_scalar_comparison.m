function plot_thf_scalar_comparison(options)
   %PLOT_THF_SCALAR_COMPARISON Three-way THF comparison: BR, scalar-exchange, MO.
   %
   %  plot_thf_scalar_comparison()
   %  plot_thf_scalar_comparison(datafile='/path/to/thf_validation_results_diagnostic.mat')
   %  plot_thf_scalar_comparison(group_idx=1:2)
   %
   % Loads results saved by run_thf_cases(output_profile='diagnostic') and
   % compares three turbulent heat flux estimates for each site/year group:
   %
   %   1. Bulk-Richardson (production scheme, single momentum roughness z0)
   %   2. Scalar exchange (BR stability factor + Andreas (2002) z0h/z0q)
   %   3. Monin-Obukhov (reference scheme)
   %
   % Produces one figure per site/year group with a 2 × 2 panel layout:
   %
   %   (1,1)  SHF 7-day moving-mean time series (all three schemes)
   %   (1,2)  LHF 7-day moving-mean time series (all three schemes)
   %   (2,1)  SHF scatter: BR vs scalar-exchange
   %   (2,2)  LHF scatter: BR vs scalar-exchange
   %
   % The scatter panels focus on BR vs. scalar-exchange because the main
   % question is whether correcting for separate scalar roughness lengths
   % materially alters the fluxes relative to the production scheme. The MO
   % scheme appears in the time-series panels for context.
   %
   % Prerequisite: run run_thf_cases(output_profile='diagnostic') first.
   %
   % See also: run_thf_cases, plot_thf_cases, validate_thf_cases

   arguments
      options.datafile (1,:) char = ''
      options.group_idx (1,:) double = []
   end

   % Default data file.
   if isempty(options.datafile)
      options.datafile = fullfile(fileparts(mfilename('fullpath')), 'data', ...
         'thf_validation_results_diagnostic.mat');
   end

   if ~isfile(options.datafile)
      error('plot_thf_scalar_comparison:missingFile', ...
         ['Data file not found:\n  %s\n' ...
         'Run run_thf_cases(output_profile=''diagnostic'') first.'], ...
         options.datafile);
   end

   loaded = load(options.datafile, 'cases', 'run_data');
   cases   = loaded.cases;
   run_data = loaded.run_data;

   figdir = fullfile(fileparts(mfilename('fullpath')), 'figures');
   if ~exist(figdir, 'dir')
      mkdir(figdir);
   end

   % Group cases by site + simulation year span.
   keys = arrayfun(@(c) ...
      sprintf('%s_%s', c.site, strjoin(string(c.simyears), '-')), cases, ...
      'UniformOutput', false);
   groups = unique(keys);

   if ~isempty(options.group_idx)
      groups = groups(options.group_idx);
      fprintf('Plotting groups: %s\n', strjoin(groups, ', '));
   end

   for g = 1:numel(groups)
      gdata = run_data(strcmp(keys, groups{g}));
      gname  = groups{g};
      gtitle = strrep(gname, '_', ' ');

      fig = plot_group(gdata, gtitle);
      fname = fullfile(figdir, [gname '_scalar_comparison.png']);
      saveas(fig, fname);
      fprintf('Saved: %s\n', fname);
   end
end

% =========================================================================
% Per-group figure
% =========================================================================

function fig = plot_group(gdata, gtitle)
   %PLOT_GROUP One figure: 2×2 panel for one site/year group.

   % Color scheme (consistent with plot_thf_cases where overlapping).
   C_br     = [0.07 0.58 0.55];   % teal  — bulk-Richardson production
   C_scalar = [0.50 0.20 0.70];   % purple — scalar-exchange experiment
   C_mo     = [0.85 0.50 0.10];   % amber  — Monin-Obukhov

   % Separate bulk_richardson and monin_obukhov runs.
   d_br = [];
   d_mo = [];
   for di = 1:numel(gdata)
      d = gdata{di};
      if strcmp(d.c.scheme, 'bulk_richardson')
         d_br = d;
      else
         d_mo = d;
      end
   end

   if isempty(d_br)
      warning('plot_thf_scalar_comparison:missingBR', ...
         'No bulk_richardson run found for group "%s" — skipping.', gtitle);
      fig = figure('Visible', 'off');
      return
   end

   % Extract daily-mean timetables.
   shf_br     = daily_mean(d_br.ice1, 'shf');
   lhf_br     = daily_mean(d_br.ice1, 'lhf');
   shf_scalar = daily_mean(d_br.ice1, 'thf_scalar_exchange_Qh');
   lhf_scalar = daily_mean(d_br.ice1, 'thf_scalar_exchange_Qe');

   fig = figure('Name', [gtitle ' — scalar comparison'], ...
      'Position', [80 80 1400 900]);
   tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
   title(tl, [gtitle ' — scalar-exchange comparison'], 'FontWeight', 'bold');

   % -----------------------------------------------------------------------
   % (1,1) SHF time series
   % -----------------------------------------------------------------------
   ax = nexttile(tl, 1);
   hold(ax, 'on');
   plot(ax, shf_br.Time,     movmean(shf_br.val,     7, 'omitnan'), ...
      '-',  'Color', C_br,     'LineWidth', 1.8, 'DisplayName', 'Bulk-Richardson');
   plot(ax, shf_scalar.Time, movmean(shf_scalar.val, 7, 'omitnan'), ...
      '--', 'Color', C_scalar, 'LineWidth', 1.5, 'DisplayName', 'Scalar exchange');
   if ~isempty(d_mo)
      shf_mo = daily_mean(d_mo.ice1, 'shf');
      plot(ax, shf_mo.Time, movmean(shf_mo.val, 7, 'omitnan'), ...
         ':', 'Color', C_mo, 'LineWidth', 1.5, 'DisplayName', 'Monin-Obukhov');
   end
   ylabel(ax, 'SHF (W m^{-2})');
   title(ax, 'Sensible Heat Flux — 7-day mean');
   legend(ax, 'Location', 'best', 'FontSize', 11, 'Box', 'off');

   % -----------------------------------------------------------------------
   % (1,2) LHF time series
   % -----------------------------------------------------------------------
   ax = nexttile(tl, 2);
   hold(ax, 'on');
   plot(ax, lhf_br.Time,     movmean(lhf_br.val,     7, 'omitnan'), ...
      '-',  'Color', C_br,     'LineWidth', 1.8, 'DisplayName', 'Bulk-Richardson');
   plot(ax, lhf_scalar.Time, movmean(lhf_scalar.val, 7, 'omitnan'), ...
      '--', 'Color', C_scalar, 'LineWidth', 1.5, 'DisplayName', 'Scalar exchange');
   if ~isempty(d_mo)
      lhf_mo = daily_mean(d_mo.ice1, 'lhf');
      plot(ax, lhf_mo.Time, movmean(lhf_mo.val, 7, 'omitnan'), ...
         ':', 'Color', C_mo, 'LineWidth', 1.5, 'DisplayName', 'Monin-Obukhov');
   end
   ylabel(ax, 'LHF (W m^{-2})');
   title(ax, 'Latent Heat Flux — 7-day mean');
   legend(ax, 'Location', 'best', 'FontSize', 11, 'Box', 'off');

   % -----------------------------------------------------------------------
   % (2,1) SHF scatter: BR vs scalar-exchange
   % -----------------------------------------------------------------------
   ax = nexttile(tl, 3);
   scatter_with_fit(ax, shf_br.val, shf_scalar.val, C_scalar, ...
      'Bulk-Richardson SHF (W m^{-2})', 'Scalar exchange SHF (W m^{-2})', ...
      'SHF: BR vs scalar-exchange');

   % -----------------------------------------------------------------------
   % (2,2) LHF scatter: BR vs scalar-exchange
   % -----------------------------------------------------------------------
   ax = nexttile(tl, 4);
   scatter_with_fit(ax, lhf_br.val, lhf_scalar.val, C_scalar, ...
      'Bulk-Richardson LHF (W m^{-2})', 'Scalar exchange LHF (W m^{-2})', ...
      'LHF: BR vs scalar-exchange');
end

% =========================================================================
% Helpers
% =========================================================================

function tt = daily_mean(ice1, varname)
   %DAILY_MEAN Return a struct with Time and val for a daily-mean variable.

   if ~ismember(varname, ice1.Properties.VariableNames)
      tt = struct('Time', NaT(0,1), 'val', zeros(0,1));
      return
   end
   d = retime(ice1(:, varname), 'daily', 'mean');
   tt = struct('Time', d.Time, 'val', d.(varname));
end

function scatter_with_fit(ax, x, y, col, xlab, ylab, ptitle)
   %SCATTER_WITH_FIT Scatter + 1:1 line + best-fit line + R².

   keep = isfinite(x) & isfinite(y);
   x = x(keep);
   y = y(keep);

   hold(ax, 'on');

   if numel(x) < 2
      title(ax, ptitle);
      xlabel(ax, xlab); ylabel(ax, ylab);
      return
   end

   scatter(ax, x, y, 14, col, 'filled', 'MarkerFaceAlpha', 0.25);

   % 1:1 reference
   lims = [min([x; y]), max([x; y])];
   plot(ax, lims, lims, 'k--', 'LineWidth', 1.0, 'DisplayName', '1:1');

   % Best-fit line
   cf = polyfit(x, y, 1);
   xf = linspace(lims(1), lims(2), 200)';
   r2 = corr(x, y) ^ 2;
   plot(ax, xf, polyval(cf, xf), '-', 'Color', col, 'LineWidth', 2.0, ...
      'DisplayName', sprintf('fit  R²=%.3f', r2));

   % Bias annotation
   bias = mean(y - x);
   text(ax, 0.05, 0.93, sprintf('bias = %+.2f W m^{-2}', bias), ...
      'Units', 'normalized', 'FontSize', 11);

   xlim(ax, lims); ylim(ax, lims);
   xlabel(ax, xlab); ylabel(ax, ylab);
   title(ax, ptitle);
   legend(ax, 'Location', 'southeast', 'FontSize', 10, 'Box', 'off');
end
