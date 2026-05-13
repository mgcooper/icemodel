function f = plotscatter(case_id, kwargs)
   %PLOTSCATTER Plot target-versus-candidate scatter panels for site cases.
   %
   %  f = icemodel.verification.plotscatter("wfj", candidate=candidate)
   %
   % This companion to plotcase is intentionally limited to time-series site
   % cases. Current synthetic-process cases such as Colbeck are better reviewed
   % as time-series/process panels until analytical profile or flux targets are
   % staged explicitly.

   arguments
      case_id (1, :) string
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = "test"
      kwargs.variables string = strings(0, 1)
      kwargs.candidate = []
      kwargs.candidate_file (1, 1) string = ""
      kwargs.visible (1, 1) string = "on"
      kwargs.output_file (1, 1) string = ""
   end

   manifest = icemodel.verification.loadmanifest(case_id, ...
      "evaluation_data_root", kwargs.evaluation_data_root, ...
      "icemodel_config_casename", kwargs.icemodel_config_casename);
   targets = icemodel.verification.helpers.loadArtifact( ...
      manifest.evaluation_path, "targets");
   if ~isfield(targets, 'format') && isfield(targets, 'numerical_summa')
      % Multi-source schema (e.g. colbeck1976): default to numerical_summa.
      targets = targets.numerical_summa;
   end
   candidate = icemodel.verification.helpers.resolveCandidateBundle(manifest, ...
      "candidate", kwargs.candidate, ...
      "candidate_file", kwargs.candidate_file);

   if targets.format ~= "timeseries"
      error('scatter plots are only supported for timeseries verification cases')
   end

   if isempty(kwargs.variables)
      variable_names = manifest.comparison_variables;
   else
      variable_names = reshape(kwargs.variables, [], 1);
   end

   n_cols = ceil(sqrt(numel(variable_names)));
   n_rows = n_cols;

   f = figure('Visible', kwargs.visible, 'Color', 'w');
   setFigureSize(f, max(720, 300 * n_cols), max(600, 300 * n_rows));

   tl = tiledlayout(f, 'flow', 'TileSpacing', 'compact', 'Padding', 'compact');
   % tl = tiledlayout(f, n_rows, n_cols, ...
   %    'TileSpacing', 'compact', 'Padding', 'compact');

   title(tl, sprintf('%s (%s) scatter', manifest.case_id, ...
      manifest.dataset_family), 'Interpreter', 'none', ...
      'FontSize', 14, 'FontWeight', 'normal')

   for n = 1:numel(variable_names)
      ax = nexttile(tl);
      plotOneVariable(ax, targets.data, candidate.data, variable_names(n));
   end
   % for n = numel(variable_names)+1:n_rows*n_cols
   %    axis(nexttile(tl), 'off')
   % end

   if ~isblanktext(kwargs.output_file)
      outdir = fileparts(kwargs.output_file);
      if ~isblanktext(outdir)
         icemodel.helpers.ensureDirExists(outdir);
      end
      exportgraphics(f, kwargs.output_file, 'Resolution', 140)
   end
end

function plotOneVariable(ax, target_tt, candidate_tt, varname)
   %PLOTONEVARIABLE Plot one finite paired scatter comparison.

   missing_target = ~ismember(varname, target_tt.Properties.VariableNames);
   missing_candidate = ~ismember(varname, candidate_tt.Properties.VariableNames);
   if missing_target || missing_candidate
      noDataPanel(ax, varname, missing_target, missing_candidate)
      return
   end

   sync_tt = synchronize(target_tt(:, varname), candidate_tt(:, varname), ...
      'intersection', 'first');
   sync_tt.Properties.VariableNames = {'target', 'candidate'};

   icemodel.plot.scatterplot(sync_tt.target, sync_tt.candidate, ...
      axes=ax, ...
      display_name="candidate", ...
      colors=[0 0.45 0.74], ...
      x_label="targets", ...
      y_label="candidate", ...
      marker_size=12, ...
      line_width=1.4, ...
      legend_location="best", ...
      legend_font_size=12);
   title(ax, strrep(varname + " scatter", '_', '\_'), ...
      'FontSize', 12, 'FontWeight', 'bold')
   formatAxes(ax);
end

function noDataPanel(ax, varname, missing_target, missing_candidate)
   %NODATAPANEL Render a visible "no data" tile and warn loudly.
   %
   % Replaces the prior silent axis-off behavior so missing variables
   % surface in both the figure and the command window.

   if missing_target && missing_candidate
      reason = 'target and candidate';
   elseif missing_target
      reason = 'target';
   else
      reason = 'candidate';
   end
   warning('icemodel:verification:plotscatter:noData', ...
      '%s scatter unavailable: %s missing', varname, reason);

   title(ax, strrep(varname + " scatter", '_', '\_'), ...
      'FontSize', 12, 'FontWeight', 'bold')
   text(ax, 0.5, 0.5, sprintf('no data\n(missing: %s)', reason), ...
      'Units', 'normalized', 'HorizontalAlignment', 'center', ...
      'VerticalAlignment', 'middle', 'FontSize', 11, 'Color', [0.4 0.4 0.4])
   ax.XTick = [];
   ax.YTick = [];
   ax.Box = 'on';
end

function setFigureSize(f, width, height)
   %SETFIGURESIZE Use stable pixel sizes for exported verification figures.

   f.Units = 'pixels';
   f.Position(3:4) = [width height];
end

function formatAxes(ax)
   %FORMATAXES Apply compact typography without making labels unreadable.

   ax.FontSize = 12;
   ax.LineWidth = 0.8;
   ax.XLabel.FontSize = 12;
   ax.YLabel.FontSize = 12;
end
