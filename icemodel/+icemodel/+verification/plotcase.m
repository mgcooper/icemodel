function f = plotcase(case_id, kwargs)
   %PLOTCASE Plot staged verification data without requiring model output.
   %
   %  f = icemodel.verification.plotcase("cdp")
   %  f = icemodel.verification.plotcase("wfj", source="reference")
   %  f = icemodel.verification.plotcase("colbeck1976", source="compare", ...
   %     candidate_file="candidate.mat", visible="off", output_file="plot.png")
   %
   % Inputs
   %  case_id                    Staged verification case id.
   %  evaluation_data_root       Base evaluation-data root. When blank, the
   %                             path is resolved from icemodel.config.
   %  icemodel_config_casename   Config casename used to resolve the default
   %                             evaluation-data root without mutating config.
   %  source                     "targets", "reference", or "compare".
   %  variables                  Optional variable subset. Defaults to the
   %                             manifest comparison variables.
   %  candidate                  Optional in-memory candidate bundle used when
   %                             source="compare".
   %  candidate_file             Optional MAT file containing `candidate` or
   %                             `reference`.
   %  visible                    Figure Visible property.
   %  output_file                Optional image path to export.
   %
   % Outputs
   %  f   Figure handle containing the staged-data plot.
   %
   % Role
   %  This is a normal verification workflow entry point for visual inspection
   %  of staged data and candidate comparisons.

   arguments
      case_id (1, :) string
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = "test"
      kwargs.source (1, 1) string ...
         {mustBeMember(kwargs.source, ["targets", "reference", "compare"])} ...
         = "targets"
      kwargs.variables string = strings(0, 1)
      kwargs.candidate = []
      kwargs.candidate_file (1, 1) string = ""
      kwargs.visible (1, 1) string = "on"
      kwargs.output_file (1, 1) string = ""
   end

   % Resolve the manifest and load the target and reference bundles. The
   % selected source below controls which bundle is plotted.
   manifest = icemodel.verification.loadmanifest(case_id, ...
      "evaluation_data_root", kwargs.evaluation_data_root, ...
      "icemodel_config_casename", kwargs.icemodel_config_casename);
   targets = icemodel.verification.helpers.loadArtifact( ...
      manifest.evaluation_path, "targets");
   reference = icemodel.verification.helpers.loadArtifact( ...
      manifest.reference_path, "reference");

   % Select the primary and optional secondary bundle. This branch is by plot
   % intent, not dataset family, so new artifact families can reuse it.
   switch kwargs.source
      case "targets"
         primary = targets;
         secondary = [];
         labels = ["targets", ""];
      case "reference"
         primary = reference;
         secondary = [];
         labels = ["reference", ""];
      case "compare"
         primary = targets;
         secondary = icemodel.verification.helpers.resolveCandidateBundle( ...
            manifest, ...
            "candidate", kwargs.candidate, ...
            "candidate_file", kwargs.candidate_file);
         labels = ["targets", "candidate"];
   end

   % Default to the manifest comparison variables so plots and metric tables use
   % the same target surface unless the caller asks for a subset.
   if isempty(kwargs.variables)
      variable_names = manifest.comparison_variables;
   else
      variable_names = reshape(kwargs.variables, [], 1);
   end

   % Dispatch by artifact format. Timeseries and experiment bundles share the
   % same high-level plotting contract but need different panel layouts.
   f = figure('Visible', kwargs.visible, 'Color', 'w');
   switch primary.format
      case "timeseries"
         plotTimeseriesBundle(f, primary.data, secondary, variable_names, ...
            labels, manifest);
      case "experiment_bundle"
         plotExperimentBundle(f, primary.experiments, secondary, ...
            variable_names, labels, manifest);
      otherwise
         close(f)
         error('unsupported verification plot format: %s', primary.format)
   end

   % Optional export keeps direct plotting useful in tests and visual review
   % without forcing every caller to manage output folders.
   if ~isblanktext(kwargs.output_file)
      outdir = fileparts(kwargs.output_file);
      if ~isblanktext(outdir)
         icemodel.helpers.ensureDirExists(outdir);
      end
      exportgraphics(f, kwargs.output_file, 'Resolution', 140)
   end
end

function plotTimeseriesBundle(f, primary_tt, secondary, variable_names, ...
      labels, manifest)
   %PLOTTIMESERIESBUNDLE Plot one timetable bundle in stacked panels.

   setFigureSize(f, 1200, max(620, 150 * numel(variable_names)));

   % One row per variable keeps site plots comparable across ESM-SnowMIP cases.
   tl = tiledlayout(f, numel(variable_names), 1, ...
      'TileSpacing', 'compact', 'Padding', 'compact');
   title(tl, sprintf('%s (%s)', manifest.case_id, manifest.dataset_family), ...
      'Interpreter', 'none', 'FontSize', 14, 'FontWeight', 'normal')

   for n = 1:numel(variable_names)
      ax = nexttile(tl);
      varname = variable_names(n);
      if ~ismember(varname, primary_tt.Properties.VariableNames)
         continue
      end

      % Plot the primary bundle first. The shared plotting helper adds markers
      % for sparse observations so isolated ESM samples remain visible.
      icemodel.plot.timeseries(primary_tt.Time, primary_tt.(varname), ...
         axes=ax, display_name=labels(1), line_style='-', ...
         color=[0 0 0], line_width=1.0);

      % Overlay candidate/reference data only when a secondary bundle exists and
      % contains the same variable.
      add_secondary = ~isempty(secondary) ...
         && ismember(varname, secondary.data.Properties.VariableNames);

      if add_secondary
         icemodel.plot.timeseries(secondary.data.Time, ...
            secondary.data.(varname), axes=ax, display_name=labels(2), ...
            line_style='--', color=[0 0.45 0.74], line_width=1.0);
      end

      ylabel(ax, strrep(varname, '_', '\_'))
      grid(ax, 'on')
      formatAxes(ax);
      if n == 1
         legend(ax, 'Location', 'best', 'FontSize', 12, 'Box', 'on')
      end
   end
end

function plotExperimentBundle(f, primary, secondary, variable_names, ...
      labels, manifest)
   %PLOTEXPERIMENTBUNDLE Plot one experiment grid such as Colbeck exp1-exp3.

   % Experiment bundles use rows for experiments and columns for variables so
   % process-case differences are visible at a glance.
   [exp_names, exp_values] = deal(fieldnames(primary), struct2cell(primary));

   setFigureSize(f, 1180, max(760, 230 * numel(exp_names)));

   % Create the tiled layout and add a title.
   tl = tiledlayout(f, numel(exp_names), numel(variable_names), ...
      'TileSpacing', 'compact', 'Padding', 'compact');
   title(tl, sprintf('%s (%s)', manifest.case_id, manifest.dataset_family), ...
      'Interpreter', 'none', 'FontSize', 14, 'FontWeight', 'normal')

   % Validate secondary experiment names before plotting to avoid silently
   % comparing the wrong experiment rows.
   if ~isempty(secondary)
      [secondary_names, secondary_values] = deal( ...
         fieldnames(secondary.experiments), struct2cell(secondary.experiments));

      if ~isequal(exp_names, secondary_names)
         error('candidate experiment names do not match target experiment names')
      end
   else
      secondary_values = {};
   end

   % Fill the experiment-by-variable grid in the same order used by comparecase.
   for n = 1:numel(exp_names)
      tt_primary = exp_values{n};

      for m = 1:numel(variable_names)
         ax = nexttile(tl);
         varname = variable_names(m);

         if ismember(varname, tt_primary.Properties.VariableNames)
            icemodel.plot.timeseries(tt_primary.Time, tt_primary.(varname), ...
               axes=ax, display_name=labels(1), line_style='-', ...
               color=[0 0 0], line_width=1.0);
         end

         if ~isempty(secondary_values)
            tt_secondary = secondary_values{n};

            if ismember(varname, tt_secondary.Properties.VariableNames)
               icemodel.plot.timeseries(tt_secondary.Time, ...
                  tt_secondary.(varname), axes=ax, display_name=labels(2), ...
                  line_style='--', color=[0 0.45 0.74], line_width=1.0);
            end
         end

         if n == 1
            title(ax, strrep(varname, '_', '\_'), ...
               'FontSize', 14, 'FontWeight', 'bold')
         end
         if m == 1
            ylabel(ax, strrep(exp_names{n}, '_', '\_'))
         end
         if n == 1 && m == 1
            legend(ax, 'Location', 'best', 'FontSize', 12, 'Box', 'on')
         end
         grid(ax, 'on')
         formatAxes(ax);
      end
   end
end

function setFigureSize(f, width, height)
   %SETFIGURESIZE Use stable pixel sizes for exported verification figures.

   f.Units = 'pixels';
   f.Position(3:4) = [width height];
end

function formatAxes(ax)
   %FORMATAXES Apply compact typography without making labels unreadable.

   ax.FontSize = 14;
   ax.LineWidth = 0.8;
   ax.XLabel.FontSize = 14;
   ax.YLabel.FontSize = 14;
   ax.Title.FontSize = 14;
   ax.Title.FontWeight = 'bold';
end
