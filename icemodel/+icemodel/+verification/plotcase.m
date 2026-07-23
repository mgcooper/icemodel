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
%  data_root                  Whole data tree containing eval/ and input/.
%  evaluation_data_root       Base evaluation-data root. When blank, the
%                             repo-local data/eval tree is used.
%  input_data_root            Optional paired input-data root for staged
%                             met/userdata artifacts.
%  icemodel_config_casename   Config casename used to resolve the default
%                             evaluation-data root without mutating config.
%  dataset_family             Optional family filter for shared case ids.
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
      kwargs.data_root (1, 1) string = ""
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.input_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = ""
      kwargs.dataset_family (1, 1) string = ""
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
      "data_root", kwargs.data_root, ...
      "evaluation_data_root", kwargs.evaluation_data_root, ...
      "input_data_root", kwargs.input_data_root, ...
      "icemodel_config_casename", kwargs.icemodel_config_casename, ...
      "dataset_family", kwargs.dataset_family);
   % The eval target is a forcing-agnostic, data-only observations.mat bundle
   % referenced via evaluation_path (ESM-SnowMIP, SUMup, freshly staged PROMICE,
   % and the Snow/Colbeck evaluation.mat). Only PROMICE fixtures staged before
   % the observations.mat contract lack that file; they fall back to
   % reconstituting the PROMICE-obs target on demand from the staged per-year
   % userdata files. The default reference is resolved from the staged model
   % sources the manifest actually declares.
   if isfield(manifest, 'evaluation_path') ...
         && strlength(string(manifest.evaluation_path)) > 0 ...
         && isfile(manifest.evaluation_path)
      targets = icemodel.verification.helpers.loadArtifact( ...
         manifest.evaluation_path, "targets");
   else
      targets = icemodel.verification.helpers.loadColocatedData( ...
         manifest, "promice");
   end
   if ~isfield(targets, 'format') && isfield(targets, 'numerical_summa')
      % Multi-source schema (e.g. colbeck1976): default to numerical_summa.
      targets = targets.numerical_summa;
   end
   reference = icemodel.verification.helpers.resolveCandidateBundle(manifest);

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
      case "subsurface_profile_bundle"
         plotProfileBundle(f, primary.data, secondary, variable_names, ...
            labels, manifest);
      case "retmip_protocol_bundle"
         plotRetmipProtocolBundle(f, primary.data, secondary, variable_names, ...
            labels, manifest);
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

      % Plot the primary bundle first using the verification-wide line-only
      % contract; source line styles distinguish overlaid bundles.
      icemodel.plot.timeseries(primary_tt.Time, primary_tt.(varname), ...
         axes=ax, display_name=labels(1), line_style='-', ...
         color=[0 0 0], line_width=1.0, marker_style='none');

      % Overlay candidate/reference data only when a secondary bundle exists and
      % contains the same variable.
      add_secondary = ~isempty(secondary) ...
         && ismember(varname, secondary.data.Properties.VariableNames);

      if add_secondary
         icemodel.plot.timeseries(secondary.data.Time, ...
            secondary.data.(varname), axes=ax, display_name=labels(2), ...
            line_style='--', color=[0 0.45 0.74], line_width=1.0, ...
            marker_style='none');
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
               color=[0 0 0], line_width=1.0, marker_style='none');
         end

         if ~isempty(secondary_values)
            tt_secondary = secondary_values{n};

            if ismember(varname, tt_secondary.Properties.VariableNames)
               icemodel.plot.timeseries(tt_secondary.Time, ...
                  tt_secondary.(varname), axes=ax, display_name=labels(2), ...
                  line_style='--', color=[0 0.45 0.74], line_width=1.0, ...
                  marker_style='none');
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
            legend(ax, 'Location', 'eastoutside', 'FontSize', 10, 'Box', 'on')
         end
         grid(ax, 'on')
         formatAxes(ax);
      end
   end
end

function plotProfileBundle(f, primary, secondary, variable_names, labels, ...
      manifest)
   %PLOTPROFILEBUNDLE Plot one firn profile/point bundle in stacked panels.

   setFigureSize(f, 1100, max(620, 170 * numel(variable_names)));
   tl = tiledlayout(f, numel(variable_names), 1, ...
      'TileSpacing', 'compact', 'Padding', 'compact');
   title(tl, sprintf('%s (%s)', manifest.case_id, manifest.dataset_family), ...
      'Interpreter', 'none', 'FontSize', 14, 'FontWeight', 'normal')

   for n = 1:numel(variable_names)
      ax = nexttile(tl);
      varname = variable_names(n);
      plotted = plotProfileWithSharedHelper( ...
         ax, primary, secondary, varname, labels);
      if ~plotted
         plotted = plotIntervalWithSharedHelper( ...
            ax, primary, secondary, varname, labels);
      end
      if ~plotted
         text(ax, 0.5, 0.5, 'No depth/time axis available', ...
            'Units', 'normalized', 'HorizontalAlignment', 'center')
         axis(ax, 'off')
      end

      grid(ax, 'on')
      formatAxes(ax);
      if n == 1
         legend(ax, 'Location', 'best', 'FontSize', 12, 'Box', 'on')
      end
   end
end

function used = plotProfileWithSharedHelper(ax, primary, secondary, varname, ...
      labels)
   %PLOTPROFILEWITHSHAREDHELPER Use the canonical depth-profile plotter.
   data = {};
   names = strings(1, 0);
   if isfield(primary, varname)
      [primary_data, primary_names] = profilePlotGroups( ...
         primary.(char(varname)), varname, labels(1));
      data = [data, primary_data];
      names = [names, primary_names];
   end
   [has_secondary, secondary_value] = secondaryValueForVariable( ...
      secondary, varname);
   if has_secondary
      [secondary_data, secondary_names] = profilePlotGroups( ...
         secondary_value, varname, labels(2));
      data = [data, secondary_data];
      names = [names, secondary_names];
   end
   if isempty(data)
      used = false;
      return
   end

   out = icemodel.plot.profile(data, varname, names=names, axes=ax);
   used = any(out.plotted);
end

function [data, names] = profilePlotGroups(value, varname, label)
   %PROFILEPLOTGROUPS Keep dated profile rows on separate plotted lines.
   data = {};
   names = strings(1, 0);
   if ~istable(value)
      profile = profileTableForVariable(value, varname);
      if ~isempty(profile)
         data = {profile};
         names = label;
      end
      return
   end

   groups = icemodel.verification.helpers.profileGroups(value);
   for n = 1:numel(groups)
      profile = profileTableForVariable(groups(n).data, varname);
      if isempty(profile)
         continue
      end
      data{end+1} = profile; %#ok<AGROW>
      if isnat(groups(n).datetime) || isscalar(groups)
         names(end+1) = label; %#ok<AGROW>
      else
         names(end+1) = label + " " ...
            + string(groups(n).datetime, 'yyyy-MM-dd'); %#ok<AGROW>
      end
   end
end

function used = plotIntervalWithSharedHelper(ax, primary, secondary, varname, ...
      labels)
   %PLOTINTERVALWITHSHAREDHELPER Plot profile-bundle interval/time variables.
   data = {};
   names = strings(1, 0);
   if isfield(primary, varname)
      primary_series = timeTableForVariable(primary.(char(varname)), varname);
      if ~isempty(primary_series)
         data{end+1} = primary_series;
         names(end+1) = labels(1);
      end
   end
   [has_secondary, secondary_value] = secondaryValueForVariable( ...
      secondary, varname);
   if has_secondary
      secondary_series = timeTableForVariable(secondary_value, varname);
      if ~isempty(secondary_series)
         data{end+1} = secondary_series;
         names(end+1) = labels(2);
      end
   end
   if isempty(data)
      used = false;
      return
   end

   out = icemodel.plot.compareTimeseries(data, varname, names=names, axes=ax);
   used = any(out.plotted);
end

function [tf, value] = secondaryValueForVariable(secondary, varname)
   %SECONDARYVALUEFORVARIABLE Return candidate data for a profile-bundle field.
   tf = false;
   value = [];
   if isempty(secondary) || ~isfield(secondary, 'data')
      return
   end
   if isstruct(secondary.data) && isfield(secondary.data, varname)
      value = secondary.data.(char(varname));
      tf = true;
      return
   end
   if (istable(secondary.data) || istimetable(secondary.data)) ...
         && ismember(varname, string(secondary.data.Properties.VariableNames))
      value = secondary.data;
      tf = true;
   end
end

function T = profileTableForVariable(value, varname)
   %PROFILETABLEFORVARIABLE Return a table with the requested value variable.
   T = [];
   if ~isTableProfile(value)
      return
   end
   names = string(value.Properties.VariableNames);
   if ismember(varname, names)
      T = value;
   elseif ismember("value", names)
      T = renamevars(value, "value", varname);
   end
end

function tf = isTableProfile(value)
   %ISTABLEPROFILE True for profile values handled by icemodel.plot.profile.
   tf = false;
   if ~(istable(value) || istimetable(value))
      return
   end
   names = string(value.Properties.VariableNames);
   tf = any(ismember(["depth_m", "depth", "midpoint", "start_depth"], names));
end

function tt = timeTableForVariable(value, varname)
   %TIMETABLEFORVARIABLE Convert time/interval tables to plotting timetables.
   tt = timetable.empty;
   if ~(istable(value) || istimetable(value))
      return
   end

   names = string(value.Properties.VariableNames);
   value_name = valueColumnName(names, varname);
   if value_name == "" || ~isnumeric(value.(char(value_name)))
      return
   end

   if istimetable(value)
      tt = timetable(value.(char(value_name)), ...
         'RowTimes', value.Properties.RowTimes, ...
         'VariableNames', cellstr(varname));
      return
   end

   if all(ismember(["start_date", "end_date"], names)) ...
         && isdatetime(value.start_date) && isdatetime(value.end_date)
      tt = intervalTimetableForVariable(value, value_name, varname);
      return
   end

   time = tableTimeAxis(value, names);
   if isempty(time)
      return
   end
   tt = timetable(value.(char(value_name)), 'RowTimes', time, ...
      'VariableNames', cellstr(varname));
end

function tt = intervalTimetableForVariable(value, value_name, varname)
   %INTERVALTIMETABLEFORVARIABLE Preserve accumulated-observation spans.
   time = reshape([value.start_date(:) value.end_date(:)].', [], 1);
   values = reshape(repmat(value.(char(value_name)).', 2, 1), [], 1);
   tt = timetable(values, 'RowTimes', time, 'VariableNames', cellstr(varname));
   tt.Properties.UserData = struct( ...
      'interval_observation', true, ...
      'interval_variables', varname);
   tt = copyVariableMetadata(tt, value, value_name, varname);
end

function tt = copyVariableMetadata(tt, value, value_name, varname)
   %COPYVARIABLEMETADATA Preserve source units/descriptions after renaming.
   old_names = string(value.Properties.VariableNames);
   new_names = string(tt.Properties.VariableNames);
   old_idx = find(old_names == value_name, 1);
   new_idx = find(new_names == varname, 1);
   if isempty(old_idx) || isempty(new_idx)
      return
   end
   if numel(value.Properties.VariableUnits) >= old_idx
      tt.Properties.VariableUnits{new_idx} = ...
         value.Properties.VariableUnits{old_idx};
   end
   if numel(value.Properties.VariableDescriptions) >= old_idx
      tt.Properties.VariableDescriptions{new_idx} = ...
         value.Properties.VariableDescriptions{old_idx};
   end
end

function name = valueColumnName(names, varname)
   %VALUECOLUMNNAME Resolve the plotted value column for a profile-bundle table.
   if ismember(varname, names)
      name = varname;
   elseif ismember("value", names)
      name = "value";
   else
      name = "";
   end
end

function time = tableTimeAxis(tbl, names)
   %TABLETIMEAXIS Resolve a representative time for point/interval rows.
   time = [];
   if all(ismember(["start_date", "end_date"], names)) ...
         && isdatetime(tbl.start_date) && isdatetime(tbl.end_date)
      time = tbl.start_date + (tbl.end_date - tbl.start_date) ./ 2;
      return
   end
   candidates = ["time", "Time", "datetime", "date"];
   hit = candidates(find(ismember(candidates, names), 1));
   if ~isempty(hit) && isdatetime(tbl.(char(hit)))
      time = tbl.(char(hit));
   end
end

function plotRetmipProtocolBundle(f, primary, secondary, variable_names, ...
      labels, manifest)
   %PLOTRETMIPPROTOCOLBUNDLE Plot RetMIP protocol surface/profile axes.

   setFigureSize(f, 1100, max(620, 170 * numel(variable_names)));
   tl = tiledlayout(f, numel(variable_names), 1, ...
      'TileSpacing', 'compact', 'Padding', 'compact');
   title(tl, sprintf('%s (%s)', manifest.case_id, manifest.dataset_family), ...
      'Interpreter', 'none', 'FontSize', 14, 'FontWeight', 'normal')

   for n = 1:numel(variable_names)
      ax = nexttile(tl);
      varname = variable_names(n);
      [x, y, xlabel_text] = retmipProtocolPlotValues(primary, varname);
      if ~isempty(y)
         plot(ax, x, y, 'k-', 'DisplayName', labels(1), 'LineWidth', 1.0);
      end

      if ~isempty(secondary) && isstruct(secondary.data)
         [x2, y2] = retmipProtocolPlotValues(secondary.data, varname);
         if ~isempty(y2)
            hold(ax, 'on')
            plot(ax, x2, y2, '--', 'Color', [0 0.45 0.74], ...
               'DisplayName', labels(2), 'LineWidth', 1.0);
            hold(ax, 'off')
         end
      end

      ylabel(ax, strrep(varname, '_', '\_'))
      xlabel(ax, xlabel_text)
      grid(ax, 'on')
      formatAxes(ax);
      if n == 1
         legend(ax, 'Location', 'best', 'FontSize', 12, 'Box', 'on')
      end
   end
end

function [x, y, xlabel_text] = retmipProtocolPlotValues(bundle, varname)
   %RETMIPPROTOCOLPLOTVALUES Extract x/y vectors from a RetMIP bundle.
   x = [];
   y = [];
   xlabel_text = "sample";
   if ~isstruct(bundle)
      return
   end
   if isfield(bundle, 'surface')
      [x, y] = retmipTablePlotValues(bundle.surface, retmipAliases(varname));
      if ~isempty(y)
         xlabel_text = "time";
         return
      end
   end
   if ~isfield(bundle, 'profiles') || ~isstruct(bundle.profiles)
      return
   end
   profiles = bundle.profiles;
   names = fieldnames(profiles);
   for k = 1:numel(names)
      [x, y] = retmipTablePlotValues(profiles.(names{k}), ...
         retmipAliases(varname));
      if ~isempty(y)
         xlabel_text = "depth";
         return
      end
   end
end

function [x, y] = retmipTablePlotValues(tbl, aliases)
   %RETMIPTABLEPLOTVALUES Extract x/y from one RetMIP table.
   x = [];
   y = [];
   if ~(istable(tbl) || istimetable(tbl))
      return
   end
   names = string(tbl.Properties.VariableNames);
   hit = aliases(find(ismember(aliases, names), 1));
   if isempty(hit) || ~isnumeric(tbl.(char(hit)))
      return
   end
   y = double(tbl.(char(hit))(:));
   if istimetable(tbl)
      x = tbl.Properties.RowTimes;
   else
      time_hit = tableTimeVariable(tbl, names);
      if strlength(time_hit) > 0
         x = tbl.(char(time_hit));
         return
      end
      depth = ["depth", "z", "depth_m"];
      depth_hit = depth(find(ismember(depth, names), 1));
      if ~isempty(depth_hit)
         x = tbl.(char(depth_hit));
      else
         x = (1:numel(y))';
      end
   end
end

function hit = tableTimeVariable(tbl, names)
   %TABLETIMEVARIABLE Return the first datetime table variable, preferring time.
   hit = "";
   preferred = ["time", "Time"];
   preferred_hit = preferred(find(ismember(preferred, names), 1));
   if ~isempty(preferred_hit) && isdatetime(tbl.(char(preferred_hit)))
      hit = preferred_hit;
      return
   end
   for n = 1:numel(names)
      if isdatetime(tbl.(char(names(n))))
         hit = names(n);
         return
      end
   end
end

function aliases = retmipAliases(varname)
   %RETMIPALIASES Return RetMIP protocol/profile column aliases.
   switch string(varname)
      case "density"
         aliases = ["density", "density_kgm3", "rho"];
      case "subsurface_temperature"
         aliases = ["subsurface_temperature", "temp", "T"];
      case "lwc"
         aliases = ["lwc", "LWC", "liquid_water_content"];
      otherwise
         aliases = string(varname);
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
