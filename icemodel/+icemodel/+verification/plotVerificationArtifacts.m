function summary = plotVerificationArtifacts(kwargs)
   %PLOTVERIFICATIONARTIFACTS Visualize staged verification artifacts.
   %
   %  summary = icemodel.verification.plotVerificationArtifacts()
   %  summary = icemodel.verification.plotVerificationArtifacts( ...
   %     output_root=fullfile(icemodel.internal.fullpath("data"), ...
   %     "preview", "firn_staging"))
   %
   % Role
   %  Visual QA helper for staged verification data. It reads manifests from an
   %  eval tree, loads observations.mat targets plus staged met/userdata legs,
   %  and writes grouped figures per case so model developers can inspect every
   %  numeric channel before using the artifacts. With overwrite=true, prior
   %  PNGs are cleared only from each selected case's figure directory.
   %
   % See also: icemodel.verification.plotFirnArtifacts,
   %  icemodel.plot.forcing, icemodel.plot.compareTimeseries

   arguments
      kwargs.dataset_family (1, :) string ...
         {icemodel.verification.validators.mustBeDatasetFamilySelection} = "all"
      kwargs.case_ids (1, :) string = strings(1, 0)
      kwargs.output_root (1, 1) string = ""
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.input_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = ""
      kwargs.figure_root (1, 1) string = ""
      kwargs.save_figs (1, 1) logical = true
      kwargs.overwrite (1, 1) logical = false
      kwargs.visible (1, 1) logical = false
      kwargs.startdate = ""
      kwargs.enddate = ""
   end

   [evaluation_data_root, input_root] = resolveRoots(kwargs);
   families = expandFamilies(kwargs.dataset_family);
   figure_root = resolveFigureRoot(kwargs.figure_root, kwargs.output_root, ...
      evaluation_data_root);
   if kwargs.save_figs
      icemodel.helpers.ensureDirExists(figure_root);
   end

   case_sets = cell(1, numel(families));
   n_cases = 0;
   for k = 1:numel(families)
      cases = familyCases(families(k), kwargs, evaluation_data_root);
      case_sets{k} = cases;
      n_cases = n_cases + numel(cases);
   end

   case_tables = cell(1, n_cases);
   n_tables = 0;
   for n = 1:numel(case_sets)
      cases = case_sets{n};
      for k = 1:numel(cases)
         case_rows = plotOneCase(cases(k), input_root, ...
            hasExplicitInputRoot(kwargs), figure_root, kwargs.save_figs, ...
            kwargs.overwrite, kwargs.visible, kwargs.startdate, kwargs.enddate);
         n_tables = n_tables + 1;
         case_tables{n_tables} = case_rows;
      end
   end

   case_tables = case_tables(1:n_tables);
   if isempty(case_tables)
      summary = emptySummary();
   else
      summary = vertcat(case_tables{:});
   end
end

function [evaluation_data_root, input_root] = resolveRoots(kwargs)
   %RESOLVEROOTS Resolve paired eval/input roots from an output root or kwargs.
   [evaluation_data_root, input_root] = ...
      icemodel.verification.setup.resolveStagingRoots( ...
      output_root=kwargs.output_root, ...
      evaluation_data_root=kwargs.evaluation_data_root, ...
      input_data_root=kwargs.input_data_root, ...
      icemodel_config_casename=kwargs.icemodel_config_casename);
end

function root = resolveFigureRoot(figure_root, output_root, evaluation_data_root)
   %RESOLVEFIGUREROOT Keep generated figures beside the staged tree by default.
   if figure_root ~= ""
      root = figure_root;
   elseif output_root ~= ""
      root = fullfile(output_root, 'figures');
   else
      root = fullfile(fileparts(evaluation_data_root), ...
         'figures', 'verification_staging');
   end
end

function families = expandFamilies(families)
   %EXPANDFAMILIES Replace "all" with the active staged family list.
   if any(families == "all")
      families = icemodel.verification.namelists.datasetfamily();
   end
   families = reshape(unique(families, 'stable'), 1, []);
end

function cases = familyCases(family, kwargs, evaluation_data_root)
   %FAMILYCASES Resolve one family and apply an optional case-id subset.
   if usesDefaultManifestDiscovery(kwargs)
      cases = icemodel.verification.listcases(dataset_family=family);
   else
      cases = icemodel.verification.listcases( ...
         dataset_family=family, ...
         evaluation_data_root=evaluation_data_root, ...
         icemodel_config_casename="");
   end
   if isempty(cases) || isempty(kwargs.case_ids)
      return
   end

   ids = string({cases.case_id});
   cases = cases(ismember(ids, kwargs.case_ids));
end

function tf = usesDefaultManifestDiscovery(kwargs)
   %USESDEFAULTMANIFESTDISCOVERY True when listcases should resolve default roots.
   tf = kwargs.output_root == "" && kwargs.evaluation_data_root == "" ...
      && kwargs.input_data_root == "" ...
      && kwargs.icemodel_config_casename == "";
end

function tf = hasExplicitInputRoot(kwargs)
   %HASEXPLICITINPUTROOT True when the caller chose the input artifact root.
   tf = kwargs.output_root ~= "" || kwargs.evaluation_data_root ~= "" ...
      || kwargs.input_data_root ~= "" ...
      || kwargs.icemodel_config_casename ~= "";
end

function rows = plotOneCase(c, input_root, explicit_input_root, figure_root, ...
      save_figs, overwrite, visible, startdate, enddate)
   %PLOTONECASE Create grouped visual-QA figures for one manifest case.
   input_root = caseInputRoot(c, input_root, explicit_input_root);
   target = loadTarget(c);
   met = stagedFiles(c, input_root, "met_files", "met");
   data = stagedFiles(c, input_root, "data_files", "userdata");
   records = timeRecords(c, target, data);
   profiles = profileRecords(target);
   groups = variableGroups(records);

   case_dir = fullfile(figure_root, char(c.dataset_family), char(c.case_id));
   if save_figs
      icemodel.helpers.ensureDirExists(case_dir);
      if overwrite
         % Remove prior case PNGs so renamed or split groups cannot linger.
         prior_figures = dir(fullfile(case_dir, '*.png'));
         for n = 1:numel(prior_figures)
            delete(fullfile(prior_figures(n).folder, prior_figures(n).name))
         end
      end
   end

   if usesDedicatedColbeckPlot(c, target)
      rows = plotColbeckComparison(c, target, met, data, input_root, ...
         case_dir, save_figs, visible);
      return
   end

   rows = cell(1, numel(groups) + 3);
   n_rows = 0;
   plotted = strings(1, 0);

   if ~isempty(met.items)
      [row, vars] = plotMetForcing(c, target, met, data, case_dir, ...
         save_figs, visible, startdate, enddate);
      n_rows = n_rows + 1;
      rows{n_rows} = row;
      plotted = unique([plotted, vars], 'stable');
   end

   for k = 1:numel(groups)
      group = groups(k);
      if isempty(group.variables)
         continue
      end
      [row, vars] = plotTimeseriesGroup(c, records, group, target, met, ...
         data, case_dir, save_figs, visible, startdate, enddate);
      n_rows = n_rows + 1;
      rows{n_rows} = row;
      plotted = unique([plotted, vars], 'stable');
   end

   if ~isempty(profiles)
      [row, vars] = plotProfiles(c, profiles, target, met, data, case_dir, ...
         save_figs, visible, startdate, enddate);
      n_rows = n_rows + 1;
      rows{n_rows} = row;
      plotted = unique([plotted, vars], 'stable');
   end

   all_numeric = allNumericVariables(records, profiles);
   % Keep diagnostic divergence, screening internals, and raw measurement
   % geometry out of the catch-all panel. They remain staged and are reported
   % as intentionally unplotted rather than obscuring the physical state.
   geometry = ismember(all_numeric, ...
      ["boom_height", "stake_height", "transducer_depth"]);
   geometry = geometry | ~cellfun('isempty', regexp(cellstr(all_numeric), ...
      '^dtice\d+$', 'once'));
   intentionally_unplotted = ["sndiv", "step_magnitude", ...
      "tice10m_source", all_numeric(geometry)];
   remaining = setdiff(all_numeric, ...
      [plotted, intentionally_unplotted], 'stable');
   if ~isempty(remaining)
      group = struct('name', "other_variables", 'title', ...
         "other numeric variables", 'display_title', ...
         "other numeric variables", 'variables', remaining, ...
         'aggregation', "mean", 'frequency', "daily");
      [row, vars] = plotTimeseriesGroup(c, records, group, target, met, ...
         data, case_dir, save_figs, visible, startdate, enddate);
      n_rows = n_rows + 1;
      rows{n_rows} = row;
      plotted = unique([plotted, vars], 'stable');
   end

   if n_rows == 0
      rows = summaryRow(c, "no_plottable_data", "", ...
         strings(1, 0), all_numeric, target, met, data);
      return
   else
      rows = rows(1:n_rows);
      missing = setdiff(all_numeric, plotted, 'stable');
      rows{n_rows} = addUnplotted(rows{n_rows}, missing);
      rows = vertcat(rows{:});
   end
end

function tf = usesDedicatedColbeckPlot(c, target)
   %USESDEDICATEDCOLBECKPLOT Identify the one staged sparse benchmark view.
   payload = target.payload;
   tf = string(c.dataset_family) == "laugh_tests" ...
      && string(c.case_id) == "colbeck1976" ...
      && isfield(payload, 'numerical_summa') ...
      && isfield(payload, 'analytical_clark2017');
end

function row = plotColbeckComparison(c, target, met, data, input_root, ...
      case_dir, save_figs, visible)
   %PLOTCOLBECKCOMPARISON Reuse plotcase's sparse experiment-grid view.
   figfile = "";
   if save_figs
      figfile = fullfile(case_dir, 'colbeck_comparison.png');
   end

   % Compare the two cached Laugh-Test targets without running a model or
   % forcing sparse experiments through the generic grouped-timeseries view.
   fig = icemodel.verification.plotcase(string(c.case_id), ...
      evaluation_data_root=fileparts(string(c.family_root)), ...
      input_data_root=string(input_root), ...
      dataset_family=string(c.dataset_family), source="compare", ...
      candidate=target.payload.analytical_clark2017, ...
      visible=string(visibleState(visible)), output_file=string(figfile));
   close(fig)

   row = summaryRow(c, "colbeck_comparison", figfile, target.variables, ...
      strings(1, 0), target, met, data);
end

function input_root = caseInputRoot(c, default_input_root, explicit_input_root)
   %CASEINPUTROOT Prefer the input root resolved with the manifest case.
   input_root = default_input_root;
   if ~explicit_input_root && isfield(c, 'input_data_root') ...
         && strlength(string(c.input_data_root)) > 0
      input_root = string(c.input_data_root);
   end
end

function [row, plotted] = plotMetForcing(c, target, met, data, case_dir, ...
      save_figs, visible, startdate, enddate)
   %PLOTMETFORCING Plot the canonical model-forcing contract in compact panels.
   records = stagedTimeRecords(met, "met");
   panels = metForcingPanels(records);
   fig = figure('Name', figureName(c, "met forcing"), ...
      'Visible', visibleState(visible), 'Color', 'w');
   setFigureSize(fig, 1180, max(620, 220 * numel(panels)));
   tl = tiledlayout(fig, numel(panels), 1, ...
      'TileSpacing', 'compact', 'Padding', 'compact');
   title(tl, sprintf('%s / %s - staged met forcing', ...
      c.dataset_family, c.case_id), 'Interpreter', 'none', ...
      'FontSize', 13, 'FontWeight', 'normal');

   plotted = strings(1, 0);
   panel_axes = gobjects(numel(panels), 1);
   for k = 1:numel(panels)
      ax = nexttile(tl);
      panel_axes(k) = ax;
      vars = plotMetPanel(ax, records, panels(k), startdate, enddate);
      plotted = unique([plotted, vars], 'stable');
   end
   shareFiniteTimeExtent(panel_axes)

   figfile = exportFigure(fig, case_dir, "met_forcing", save_figs);
   close(fig)

   row = summaryRow(c, "met_forcing", figfile, plotted, strings(1, 0), ...
      target, met, data);
end

function panels = metForcingPanels(records)
   %METFORCINGPANELS Return at most five full-width canonical forcing panels.
   prototype = struct('title', "", 'left', strings(1, 0), ...
      'right', strings(1, 0), 'aggregation', "mean", ...
      'left_label', "", 'right_label', "");
   candidates = repmat(prototype, 1, 5);
   candidates(1) = struct('title', "downwelling radiation", ...
      'left', ["swd", "lwd"], 'right', strings(1, 0), ...
      'aggregation', "mean", ...
      'left_label', "downwelling radiation [W m-2]", 'right_label', "");
   candidates(2) = struct('title', "air temperature and humidity", ...
      'left', "tair", 'right', "rh", 'aggregation', "mean", ...
      'left_label', "tair [degC]", 'right_label', "rh [%]");
   candidates(3) = struct('title', "wind speed and surface pressure", ...
      'left', "wspd", 'right', "psfc", 'aggregation', "mean", ...
      'left_label', "wspd [m s-1]", 'right_label', "psfc [Pa]");
    candidates(4) = struct('title', "surface albedo", 'left', "albedo", ...
       'right', strings(1, 0), ...
       'aggregation', "shortwave_weighted_mean", ...
       'left_label', "albedo [-]", 'right_label', "");
   candidates(5) = struct('title', "precipitation", 'left', "ppt", ...
      'right', strings(1, 0), 'aggregation', "daily_total", ...
      'left_label', "ppt [m day-1]", 'right_label', "");

   present = strings(1, 0);
   if ~isempty(records)
      present = unique([records.variables], 'stable');
   end
   keep = arrayfun(@(p) any(ismember([p.left, p.right], present)), candidates);
   panels = candidates(keep);
end

function plotted = plotMetPanel(ax, records, panel, startdate, enddate)
   %PLOTMETPANEL Render one canonical forcing panel, with a paired axis if used.
   % Keep the panel identity available without repeating its y-axis label as a
   % visible title in every tile of the stacked forcing figure.
   ax.Tag = char(panel.title);
   paired = ~isempty(panel.right);
   if paired
      % Create paired rulers only when the panel actually uses its right axis.
      yyaxis(ax, 'left')
   end
   % Use line patterns consistently within a source: the first left-axis
   % scalar is solid, its paired or secondary scalar is dotted.
   left_styles = {'-', ':', '-.', '--'};
   right_styles = {':', '-.', '--', '-'};
   plotted = plotPanelVariables(ax, records, panel.left, panel.aggregation, ...
      paired, left_styles, startdate, enddate);
   if panel.left_label ~= ""
      ylabel(ax, panel.left_label, 'Interpreter', 'none')
   end

   if paired
      yyaxis(ax, 'right')
      right = plotPanelVariables(ax, records, panel.right, panel.aggregation, ...
         paired, right_styles, startdate, enddate);
      plotted = unique([plotted, right], 'stable');
      ylabel(ax, panel.right_label, 'Interpreter', 'none')
   end
   styleLegend(ax);
end

function plotted = plotPanelVariables(ax, records, variables, aggregation, ...
      qualify_variable, styles, startdate, enddate)
   %PLOTPANELVARIABLES Overlay one or more variables with source-keyed colors.
   plotted = strings(1, numel(variables));
   n_plotted = 0;
   for k = 1:numel(variables)
      varname = variables(k);
      [payloads, names, ~, albedo_sources] = ...
         recordsForVariable(records, varname);
      if isempty(payloads)
         continue
      end
      if qualify_variable || numel(variables) > 1
         % Paired scalar calls share one legend, so the variable suffix is
         % required to distinguish otherwise identical source labels.
         names = names + " - " + varname;
      end
      out = icemodel.plot.compareTimeseries(payloads, varname, axes=ax, ...
         names=names, frequency="daily", aggregation=aggregation, ...
         line_style=styles{min(k, numel(styles))}, ...
         show_legend=false, albedo_source=albedo_sources, ...
         startdate=startdate, enddate=enddate);
      if any(out.plotted)
         n_plotted = n_plotted + 1;
         plotted(n_plotted) = varname;
      end
   end
   plotted = plotted(1:n_plotted);
end

function groups = variableGroups(records)
   %VARIABLEGROUPS Return concrete plotting groups from present variables.
   if isempty(records)
      vars = strings(1, 0);
   else
      vars = unique([records.variables], 'stable');
   end
   groups = repmat(emptyGroup(), 1, 13);
   % Observed radiation has legitimate partial-hour support (especially SWD/SWU
   % at night). Plot native samples with explicit gap breaks instead of either
   % discarding partial days or biasing their means with omit-missing reduction.
   groups(1) = makeGroup("radiation_fluxes", "radiation fluxes", vars, ...
      ["swd", "swu", "lwd", "lwu"], "mean", "native");
   groups(2) = makeGroup("energy_balance", "energy-balance terms", vars, ...
      ["swn", "lwn", "thf", "netr"], "mean", "daily");
   groups(3) = makeGroup("turbulent_fluxes", "turbulent fluxes", vars, ...
      ["shf", "lhf"], "mean", "daily");
   groups(4) = makeGroup("temperature_state", "temperature state", vars, ...
      ["tair", "tsfc", "tice10m", "surface_temp_C"], "mean", "daily", ...
      pattern="^soil_temp_\d+_C$");
   groups(5) = makeGroup("mass_balance_terms", ...
      "canonical surface-mass-balance terms", vars, ...
      ["ppt", "runoff", "evap", "subl", "smb"], ...
      "daily_total", "daily");
   % SUMup SMB is a total over each source-specific observation interval, not a
   % daily increment. Keep it on a separate quantitative panel from RCM SMB.
   groups(6) = makeGroup("interval_surface_mass_balance", ...
      "observed interval surface mass balance", vars, "smb_interval", ...
      "mean", "native");
   groups(7) = makeGroup("mass_balance_components", ...
      "precipitation and meltwater components", vars, ...
      ["rainf", "snowf", "snowf_subl", "melt", "refreeze", "meltin"], ...
      "daily_total", "daily");
   groups(8) = makeGroup("mar_combined_mass_diagnostics", ...
      "MAR combined daily mass diagnostics", vars, ...
      ["subl_evap", "refreeze_deposition"], ...
      "daily_total", "daily");
   % Plot corrected physical surface state without the step-screening internals.
   % Sparse evaluation channels use available finite samples; continuous forcing
   % retains strict complete-day means in plotTimeseriesGroup.
   groups(9) = makeGroup("surface_height_depth", ...
      "surface height change, snow depth, and stores", vars, ...
      ["ablation", "surface_height", "snow_depth", "swe", ...
       "snow_depth_m", "swe_kg_m2"], "mean", "daily", ...
       pattern="^snd_.*_m$");
   groups(9).title = "surface height change, snow depth, and stores " ...
      + "(daily means; complete forcing days, available observation samples)";
   % Albedo is structurally absent without usable shortwave radiation. Weight
   % its daily mean by measured SWD while still requiring a complete native
   % timestamp grid; every other dense variable keeps strict completeness.
   groups(10) = makeGroup("surface_albedo", "surface albedo", vars, ...
      "albedo", "shortwave_weighted_mean", "daily");
   if any(ismember(["albedo", "modis"], vars))
      % MODIS is drawn as the observational albedo overlay, not its own panel.
      groups(10).variables = "albedo";
   end
   groups(11) = makeGroup("meteorological_diagnostics", ...
      "meteorological diagnostics", vars, ...
      ["rh", "wspd", "wdir", "psfc", "cfrac"], ...
      "mean", "daily");
   groups(12) = makeGroup("quality_flags", "quality flags", vars, ...
      ["surface_height_flag", "station_transition_flag", ...
      "step_detected_flag", "step_correctable_flag", ...
      "tice10m_qc_flag"], "mean", "native");
   groups(13) = makeGroup("subsurface_temperature_string", ...
      "thermistor string", vars, strings(1, 0), "mean", "native", ...
      pattern="^tice\d+$");
   groups = splitOversizedGroups(groups, 7);
end

function groups = splitOversizedGroups(groups, max_panels)
   %SPLITOVERSIZEDGROUPS Keep every full-width figure within its panel limit.
   expanded = repmat(emptyGroup(), 1, 0);
   for k = 1:numel(groups)
      group = groups(k);
      if group.name == "subsurface_temperature_string"
         % The established PROMICE QA view overlays the complete thermistor
         % string on one axes; it is not a stack of per-sensor panels.
         expanded(end + 1) = group; %#ok<AGROW>
         continue
      end
      n_parts = max(1, ceil(numel(group.variables) / max_panels));
      for part_index = 1:n_parts
         first = (part_index - 1) * max_panels + 1;
         last = min(part_index * max_panels, numel(group.variables));
         part = group;
         part.variables = group.variables(first:last);
         if n_parts > 1
            % Preserve the established first-group name for report compatibility;
            % later chunks receive deterministic suffixes and readable titles.
            part.title = group.title + " (" + part_index + " of " + n_parts + ")";
            part.display_title = group.display_title ...
               + " (" + part_index + " of " + n_parts + ")";
            if part_index > 1
               part.name = group.name + "_" + part_index;
            end
         end
         expanded(end + 1) = part; %#ok<AGROW>
      end
   end
   groups = expanded;
end

function group = emptyGroup()
   %EMPTYGROUP Prototype for a variable plotting group.
   group = struct('name', "", 'title', "", 'display_title', "", ...
      'variables', strings(1, 0), 'aggregation', "mean", ...
      'frequency', "daily");
end

function group = makeGroup(name, title_text, vars, preferred, aggregation, ...
      frequency, kwargs)
   %MAKEGROUP Intersect preferred variables and optional regex pattern.
   arguments
      name (1, 1) string
      title_text (1, 1) string
      vars (1, :) string
      preferred (1, :) string
      aggregation (1, 1) string
      frequency (1, 1) string
      kwargs.pattern (1, 1) string = ""
   end

   present = preferred(ismember(preferred, vars));
   if kwargs.pattern ~= ""
      match = vars(~cellfun('isempty', regexp(cellstr(vars), ...
         char(kwargs.pattern), 'once')));
      present = unique([present, match], 'stable');
   end
   display_title = title_text;
   title_text = readinessTitle(title_text, aggregation, frequency);
   group = struct('name', name, 'title', title_text, ...
      'display_title', display_title, 'variables', present, ...
      'aggregation', aggregation, 'frequency', frequency);
end

function title_text = readinessTitle(title_text, aggregation, frequency)
   %READINESSTITLE Distinguish plot reductions from native source support.
   if frequency == "native"
      title_text = title_text + " (native support)";
   elseif frequency ~= "daily"
      return
   elseif aggregation == "shortwave_weighted_mean"
      % Albedo reduction depends on provenance: radiometers use daily energy
      % ratios with a support gate, model states use finite-sample means, and
      % model flux ratios use daily energy ratios without that observation gate.
      title_text = title_text ...
         + " (daily source-aware means; source-specific support rules)";
   elseif aggregation == "daily_total"
      title_text = title_text + " (daily totals; complete days only)";
   else
      title_text = title_text + " (daily means; complete days only)";
   end
end

function [row, plotted] = plotTimeseriesGroup(c, records, group, target, ...
      met, data, case_dir, save_figs, visible, startdate, enddate)
   %PLOTTIMESERIESGROUP Plot one variable group as one panel per variable.
   if group.name == "subsurface_temperature_string"
      [row, plotted] = plotThermistorStringGroup(c, records, group, target, ...
         met, data, case_dir, save_figs, visible, startdate, enddate);
      return
   end

   fig = figure('Name', figureName(c, group.title), ...
      'Visible', visibleState(visible), 'Color', 'w');
   setFigureSize(fig, 1180, max(620, 220 * numel(group.variables)));
   tl = tiledlayout(fig, numel(group.variables), 1, ...
      'TileSpacing', 'compact', 'Padding', 'compact');
   title(tl, sprintf('%s / %s - %s', c.dataset_family, c.case_id, ...
      group.display_title), 'Interpreter', 'none', 'FontSize', 13, ...
      'FontWeight', 'normal');

   plotted = strings(1, numel(group.variables));
   panel_axes = gobjects(numel(group.variables), 1);
   n_plotted = 0;
   for k = 1:numel(group.variables)
      varname = group.variables(k);
      ax = nexttile(tl);
      % A nonvisual tag keeps semantic panel selection stable after removing
      % titles that merely repeated the y-axis variable name.
      ax.Tag = char(varname);
      panel_axes(k) = ax;
       [payloads, names, observations, albedo_sources] = ...
          recordsForVariable(records, varname);
      if isempty(payloads)
         noDataPanel(ax, varname, "not present");
         continue
      end
      aggregations = repmat(group.aggregation, 1, numel(payloads));
       if group.name == "surface_height_depth"
          aggregations(observations) = "mean_omitmissing";
       end
       display_unit = groupDisplayUnit(group.name, varname);
       out = icemodel.plot.compareTimeseries(payloads, varname, axes=ax, ...
          names=names, frequency=group.frequency, ...
          aggregation=aggregations, startdate=startdate, ...
          enddate=enddate, display_unit=display_unit, ...
          albedo_source=albedo_sources, show_legend=false);
       explanatory_title = variableDisplayTitle(varname);
       if explanatory_title ~= ""
          % Retain only titles that explain semantics not present in the axis
          % label; reduction/support details remain in the figure/report text.
          title(ax, explanatory_title, 'FontWeight', 'bold', ...
             'Interpreter', 'none')
       end
       if varname == "smb_interval"
          % The display alias names the source support, while this shorter axis
          % label keeps the physical quantity and unit readable in the report.
          ylabel(ax, "interval SMB [" + out.unit + "]", ...
             'Interpreter', 'none')
      end
      if any(out.plotted)
         styleLegend(ax);
         n_plotted = n_plotted + 1;
         plotted(n_plotted) = varname;
      else
         noDataPanel(ax, varname, "not numeric or all missing");
      end
   end
   shareFiniteTimeExtent(panel_axes)

   plotted = plotted(1:n_plotted);
   if group.name == "surface_albedo" ...
         && any(arrayfun(@(r) ismember("modis", r.variables), records))
      plotted = unique([plotted, "modis"], 'stable');
   end
   figfile = exportFigure(fig, case_dir, group.name, save_figs);
   close(fig)
   row = summaryRow(c, group.name, figfile, plotted, strings(1, 0), ...
      target, met, data);
end

function unit = groupDisplayUnit(group_name, varname)
   %GROUPDISPLAYUNIT Return presentation units without changing staged metadata.

   mass_groups = ["mass_balance_terms", "mass_balance_components", ...
      "mar_combined_mass_diagnostics"];
   if ismember(group_name, mass_groups) && varname ~= "smb"
      % These canonical surface-mass-balance channels are water-equivalent
      % increments; millimetres keep daily values legible across sources.
      unit = "mmWE";
   else
      unit = "";
   end
end

function label = variableDisplayTitle(varname)
   %VARIABLEDISPLAYTITLE Explain variables whose compact names hide semantics.

   switch varname
      case "snowf_subl"
         label = "net snow accumulation (+snowfall/deposition; -sublimation)";
      case "surface_height"
         label = "surface height change (+up; not snow depth)";
      case "smb_interval"
         label = "surface mass balance accumulated over observation interval";
      case "snow_depth"
         label = "";
      case "tice10m"
         label = "10 m ice/firn temperature (canonical tice10m)";
      otherwise
         label = "";
   end
end

function [row, plotted] = plotThermistorStringGroup(c, records, group, ...
      target, met, data, case_dir, save_figs, visible, startdate, enddate)
   %PLOTTHERMISTORSTRINGGROUP Reuse the accepted one-axes PROMICE QA design.
   fig = figure('Name', figureName(c, group.title), ...
      'Visible', visibleState(visible), 'Color', 'w');
   setFigureSize(fig, 1180, 620);
   tl = tiledlayout(fig, 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
   title(tl, sprintf('%s / %s - %s', c.dataset_family, c.case_id, ...
      group.display_title), 'Interpreter', 'none', 'FontSize', 13, ...
      'FontWeight', 'normal');
   ax = nexttile(tl);
   hold(ax, 'on')

   % Match the established PROMICE QA figure: every depth-tagged thermistor is
   % a thin grey diagnostic line on the same axes.
   plotted = strings(1, numel(group.variables) + 1);
   n_plotted = 0;
   for varname = reshape(group.variables, 1, [])
      [payloads, names] = recordsForVariable(records, varname);
      if isempty(payloads)
         continue
      end
      out = icemodel.plot.compareTimeseries(payloads, varname, axes=ax, ...
         names=names, frequency=group.frequency, ...
         aggregation=group.aggregation, startdate=startdate, ...
         enddate=enddate, show_legend=false);
      if ~any(out.plotted)
         continue
      end
      labels = thermistorLabels(names(out.plotted), varname);
      for n = 1:numel(out.lines)
         set(out.lines(n), 'Color', [0.6 0.6 0.6], 'LineWidth', 0.6, ...
            'DisplayName', char(labels(n)), 'Marker', 'none')
      end
      n_plotted = n_plotted + 1;
      plotted(n_plotted) = varname;
   end

   % Draw the primary QC-screened 10 m channel last and thick so it remains
   % unmistakable above the grey source string, exactly as in PROMICE QA/QC.
   [payloads, names] = recordsForVariable(records, "tice10m");
   if ~isempty(payloads)
      out = icemodel.plot.compareTimeseries(payloads, "tice10m", axes=ax, ...
         names=names, frequency=group.frequency, ...
         aggregation=group.aggregation, startdate=startdate, ...
         enddate=enddate, show_legend=false);
      if any(out.plotted)
         labels = thermistorLabels(names(out.plotted), "tice10m (PRIMARY)");
         for n = 1:numel(out.lines)
            set(out.lines(n), 'Color', [0 0 0], 'LineWidth', 2.0, ...
               'DisplayName', char(labels(n)), 'Marker', 'none')
         end
         n_plotted = n_plotted + 1;
         plotted(n_plotted) = "tice10m";
      end
   end

   % The zero guide is diagnostic only and must not add a legend entry.
   yline(ax, 0, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
   hold(ax, 'off')
   ylabel(ax, 'ice temperature [degC]')
   xlabel(ax, 'time (UTC)')
   title(ax, ['subsurface temperature: tice10m (primary) + ' ...
      'depth-tagged string (shallow->deep)'])
   if n_plotted > 0
      % This specialized legend is moved outside immediately below, so it must
      % not reserve in-axes headroom intended for the shared stacked panels.
      styleLegend(ax, false);
      lgd = findobj(fig, 'Type', 'Legend');
      lgd.Location = 'eastoutside';
      lgd.NumColumns = 1;
   else
      noDataPanel(ax, "thermistor string", "not numeric or all missing");
   end

   plotted = plotted(1:n_plotted);
   figfile = exportFigure(fig, case_dir, group.name, save_figs);
   close(fig)
   row = summaryRow(c, group.name, figfile, plotted, strings(1, 0), ...
      target, met, data);
end

function labels = thermistorLabels(source_names, channel_name)
   %THERMISTORLABELS Keep concise channel labels unless sources are duplicated.
   labels = repmat(channel_name, size(source_names));
   if numel(source_names) > 1
      labels = source_names + " - " + channel_name;
   end
end

function [row, plotted] = plotProfiles(c, profiles, target, met, data, ...
      case_dir, save_figs, visible, startdate, enddate)
   %PLOTPROFILES Plot all depth-profile variables in one figure.
   value_vars = unique([profiles.value_variables], 'stable');
   fig = figure('Name', figureName(c, "profiles"), ...
      'Visible', visibleState(visible), 'Color', 'w');
   setFigureSize(fig, 1050, max(620, 260 * numel(value_vars)));
   tl = tiledlayout(fig, numel(value_vars), 1, ...
      'TileSpacing', 'compact', 'Padding', 'compact');
   title(tl, sprintf('%s / %s - profiles', c.dataset_family, c.case_id), ...
      'Interpreter', 'none', 'FontSize', 13, 'FontWeight', 'normal')

   plotted = strings(1, numel(value_vars));
   n_plotted = 0;
   for k = 1:numel(value_vars)
      varname = value_vars(k);
      ax = nexttile(tl);
      [tables, names] = profilesForVariable(profiles, varname, c.dataset_family);
      out = icemodel.plot.profile(tables, varname, axes=ax, names=names, ...
         startdate=startdate, enddate=enddate);
      title(ax, strrep(varname, '_', '\_'), 'FontWeight', 'bold')
      if any(out.plotted)
         styleLegend(ax);
         n_plotted = n_plotted + 1;
         plotted(n_plotted) = varname;
      else
         noDataPanel(ax, varname, "profile variable not plottable");
      end
   end

   plotted = plotted(1:n_plotted);
   figfile = exportFigure(fig, case_dir, "profiles", save_figs);
   close(fig)
   row = summaryRow(c, "profiles", figfile, plotted, strings(1, 0), ...
      target, met, data);
end

function target = loadTarget(c)
   %LOADTARGET Load a case's observations/evaluation bundle if it exists.
   target = struct('format', "missing", 'payload', struct(), ...
      'variables', strings(0, 1), 'note', "missing observations.mat");
   path = string(icemodel.verification.helpers.fieldOr(c, ...
      'evaluation_path', ""));
   if path == "" || ~isfile(path)
      return
   end

   loaded = load(path, 'targets');
   if ~isfield(loaded, 'targets')
      target.note = "file has no targets variable";
      return
   end

   payload = icemodel.verification.setup.stampArtifactMetadata( ...
      loaded.targets);
   target.payload = payload;
   if isfield(payload, 'format')
      target.format = string(payload.format);
   else
      target.format = "unknown";
   end
   target.variables = targetVariables(payload);
   target.note = "";
end

function variables = targetVariables(payload)
   %TARGETVARIABLES Collect variable names from a target bundle.
   variables = tableVariables(payload);
end

function variables = tableVariables(value)
   %TABLEVARIABLES Recursively collect table/timetable variable names.
   variables = strings(1, 0);
   if istable(value) || istimetable(value)
      variables = string(value.Properties.VariableNames);
   elseif isstruct(value)
      fields = fieldnames(value);
      for k = 1:numel(fields)
         variables = unique([variables, tableVariables(value.(fields{k}))], ...
            'stable');
      end
   end
end

function files = stagedFiles(c, input_root, manifest_field, input_subdir)
   %STAGEDFILES Load staged met or userdata files declared in colocation legs.
   files = struct('sources', strings(0, 1), 'items', {{}});
   if isfield(c, 'colocation') && ~isempty(c.colocation)
      fields = string(fieldnames(c.colocation));
   else
      fields = strings(1, 0);
   end

   items = cell(1, stagedFileCount(c, manifest_field));
   n_items = 0;
   for f = reshape(fields, 1, [])
      leg = c.colocation.(f);
      if ~isstruct(leg) || ~isfield(leg, manifest_field)
         continue
      end
      rel = string(leg.(manifest_field));
      rel = rel(strlength(rel) > 0);
      for r = reshape(rel, 1, [])
         pathname = fullfile(input_root, input_subdir, r);
         if ~isfile(pathname)
            continue
         end
         n_items = n_items + 1;
         items{n_items} = struct('source', legSourceLabel(f, leg), ...
            'path', pathname, 'payload', loadStagedFile(pathname, ...
            input_subdir));
      end
   end

   if n_items == 0
      items = fallbackStagedFiles(c, input_root, manifest_field, input_subdir);
      n_items = numel(items);
   end

   files.items = items(1:n_items);
   if n_items > 0
      sources = strings(1, n_items);
      for k = 1:n_items
         sources(k) = files.items{k}.source;
      end
      files.sources = unique(sources, 'stable');
   end
end

function label = legSourceLabel(source, leg)
   %LEGSOURCELABEL Return the public manifest source id for one staged leg.
   source = string(source);
   if isstruct(leg) && isfield(leg, 'source_id') ...
         && strlength(string(leg.source_id)) > 0
      label = string(leg.source_id);
   elseif ismember(source, icemodel.verification.namelists.rcmsources())
      label = icemodel.verification.namelists.rcmProductIds(source);
   else
      label = source;
   end
end

function items = fallbackStagedFiles(c, input_root, manifest_field, input_subdir)
   %FALLBACKSTAGEDFILES Resolve atomic ESM met through the standard runtime path.
   items = {};
   if input_subdir ~= "met" || manifest_field ~= "met_files" ...
         || string(c.dataset_family) ~= "esm_snowmip"
      return
   end

   % The atomic manifest has no met_files field. Reuse the actual model option
   % resolver so nested/flat precedence, cadence, period, and enclosing-window
   % selection cannot drift from runtime behavior or admit wildcard decoys.
   paths = icemodel.verification.helpers.esmRuntimeMetFiles(c, input_root);
   paths = paths(isfile(paths));
   paths = unique(paths, 'stable');
   items = cell(1, numel(paths));
   for k = 1:numel(paths)
      items{k} = struct('source', "esm_snowmip", 'path', paths(k), ...
         'payload', loadStagedFile(paths(k), input_subdir));
   end
end

function payload = loadStagedFile(pathname, input_subdir)
   %LOADSTAGEDFILE Load the timetable payload from a staged MAT artifact.
   loaded = load(pathname);
   payload = timetable.empty;
   if input_subdir == "met" && isfield(loaded, 'met')
      payload = icemodel.verification.setup.stampArtifactMetadata(loaded.met);
   elseif input_subdir == "userdata" && isfield(loaded, 'Data')
      payload = icemodel.verification.setup.stampArtifactMetadata(loaded.Data);
   end
end

function n_files = stagedFileCount(c, manifest_field)
   %STAGEDFILECOUNT Count declared staged-file paths before loading them.
   n_files = 0;
   if ~isfield(c, 'colocation') || isempty(c.colocation)
      return
   end

   fields = string(fieldnames(c.colocation));
   for f = reshape(fields, 1, [])
      leg = c.colocation.(f);
      if ~isstruct(leg) || ~isfield(leg, manifest_field)
         continue
      end
      rel = string(leg.(manifest_field));
      n_files = n_files + nnz(strlength(rel) > 0);
   end
end

function records = timeRecords(c, target, data)
   %TIMERECORDS Build comparison records without resampled met duplication.
   records = emptyRecord();
   records(1) = [];
   records = [records, targetTimeRecords(target, string(c.dataset_family))];
   records = [records, stagedTimeRecords(data, "userdata")];
   records = addDerivedPlotVariables(records);
end

function records = addDerivedPlotVariables(records)
   %ADDDERIVEDPLOTVARIABLES Add display aliases/sums without mutating artifacts.
   for k = 1:numel(records)
      T = records(k).payload;
      names = string(T.Properties.VariableNames);

      % MAR userdata intentionally stages native radiation components only;
      % derive the model-ready balances for display with the same surface
      % emissivity convention used by icemodel.processmet. The staged file is
      % unchanged, and the in-memory metadata records that provenance.
      mar_components = ["swd", "albedo", "tsfc", "lwd"];
      mar_derived = strings(1, 3);
      num_mar_derived = 0;
      if records(k).source == "mar3.11" ...
            && all(ismember(mar_components, names))
         if ~ismember("swn", names)
            T.swn = T.swd .* (1 - T.albedo);
            num_mar_derived = num_mar_derived + 1;
            mar_derived(num_mar_derived) = "swn";
         end
         if ~ismember("lwn", names)
            T.lwn = icemodel.surface.net_longwave_radiation(T.tsfc, T.lwd);
            num_mar_derived = num_mar_derived + 1;
            mar_derived(num_mar_derived) = "lwn";
         end
         names = string(T.Properties.VariableNames);
         if ~ismember("netr", names) && all(ismember(["swn", "lwn"], names))
            T.netr = T.swn + T.lwn;
            num_mar_derived = num_mar_derived + 1;
            mar_derived(num_mar_derived) = "netr";
         end
         mar_derived = mar_derived(1:num_mar_derived);
         for name = mar_derived
            index = find(string(T.Properties.VariableNames) == name, 1);
            units = icemodel.forcing.helpers.variableUnits(name);
            T.Properties.VariableUnits{index} = units{1};
         end
         if ~isempty(mar_derived)
            metadata = T.Properties.UserData;
            metadata.display_derived_variables = mar_derived;
            metadata.display_derived_basis = struct( ...
               'swn', "swd .* (1 - albedo)", ...
               'lwn', "icemodel.surface.net_longwave_radiation(tsfc, lwd)", ...
               'netr', "swn + lwn", ...
               'scope', "display only; staged payload unchanged");
            T.Properties.UserData = metadata;
         end
      end

      names = string(T.Properties.VariableNames);
      if ~ismember("netr", names) && all(ismember(["swn", "lwn"], names))
         T.netr = T.swn + T.lwn;
         units = icemodel.forcing.helpers.variableUnits("netr");
         T.Properties.VariableUnits{end} = units{1};
      end
      names = string(T.Properties.VariableNames);
      if ~ismember("thf", names) && all(ismember(["shf", "lhf"], names))
         T.thf = T.shf + T.lhf;
         units = icemodel.forcing.helpers.variableUnits("thf");
         T.Properties.VariableUnits{end} = units{1};
      end

      % Interval-observation SMB is an accumulated source value. Give it a
      % display-only alias so the renderer cannot overlay it with daily RCM SMB.
      metadata = T.Properties.UserData;
      interval_smb = isstruct(metadata) ...
         && isfield(metadata, 'interval_observation') ...
         && logical(metadata.interval_observation) ...
         && isfield(metadata, 'interval_variables') ...
         && any(string(metadata.interval_variables) == "smb") ...
         && ismember("smb", string(T.Properties.VariableNames));
      if interval_smb
         source_index = find(string(T.Properties.VariableNames) == "smb", 1);
         T.smb_interval = T.smb;
         alias_index = width(T);
         if numel(T.Properties.VariableUnits) >= source_index
            T.Properties.VariableUnits{alias_index} = ...
               T.Properties.VariableUnits{source_index};
         end
         if numel(T.Properties.VariableDescriptions) >= source_index
            T.Properties.VariableDescriptions{alias_index} = ...
               T.Properties.VariableDescriptions{source_index};
         end
         metadata.interval_variables = unique( ...
            [string(metadata.interval_variables(:)); "smb_interval"], ...
            'stable');
         T.Properties.UserData = metadata;
      end

      % MAR calls physical snow depth SNOWD while PROMICE uses SNOW_DEPTH.
      % Add a display-only canonical alias and retain SNOWD in the in-memory
      % payload so the staged artifact name and source provenance remain intact.
      names = string(T.Properties.VariableNames);
      if ~ismember("snow_depth", names) && ismember("snowd", names)
         source_index = find(names == "snowd", 1);
         T.snow_depth = T.snowd;
         alias_index = width(T);
         if numel(T.Properties.VariableUnits) >= source_index
            T.Properties.VariableUnits{alias_index} = ...
               T.Properties.VariableUnits{source_index};
         end
         if numel(T.Properties.VariableDescriptions) >= source_index
            T.Properties.VariableDescriptions{alias_index} = ...
               T.Properties.VariableDescriptions{source_index};
         end
      end
      records(k).payload = T;
      variables = numericVariables(T);
      hidden = ["snowd", ""];
      if interval_smb
         % Keep the staged canonical field in memory but expose only its interval
         % alias to the group registry; model records still expose daily SMB.
         hidden(2) = "smb";
      end
      records(k).variables = variables( ...
         ~ismember(variables, hidden(hidden ~= "")));
   end
end

function records = targetTimeRecords(target, dataset_family)
   %TARGETTIMERECORDS Convert target time-series tables to records.
   records = emptyRecord();
   records(1) = [];
   payload = target.payload;
   if isempty(fieldnames(payload))
      return
   end
   records = recordsFromValue(payload, dataset_family, "observations");
end

function records = stagedTimeRecords(files, kind)
   %STAGEDTIMERECORDS Convert loaded met/userdata files to records.
   record_cells = cell(1, numel(files.items));
   for k = 1:numel(files.items)
      item = files.items{k};
      record_cells{k} = recordsFromValue(item.payload, ...
         string(item.source), kind);
   end
   records = concatRecords(record_cells);
end

function records = recordsFromValue(value, source, kind)
   %RECORDSFROMVALUE Convert a table-like payload to timetable records.
   records = emptyRecord();
   records(1) = [];
   if istimetable(value)
      records = makeRecord(source, kind, value);
   elseif istable(value)
      T = tableToTimeRecord(value);
      if istimetable(T)
         records = makeRecord(source, kind, T);
      end
   elseif isstruct(value)
      fields = fieldnames(value);
      record_cells = cell(1, numel(fields));
      for k = 1:numel(fields)
         record_cells{k} = recordsFromValue(value.(fields{k}), ...
            source, kind);
      end
      records = concatRecords(record_cells);
   end
end

function records = concatRecords(record_cells)
   %CONCATRECORDS Concatenate optional record arrays from recursive readers.
   records = emptyRecord();
   records(1) = [];
   if isempty(record_cells)
      return
   end
   keep = ~cellfun(@isempty, record_cells);
   record_cells = record_cells(keep);
   if ~isempty(record_cells)
      records = [record_cells{:}];
   end
end

function record = emptyRecord()
   %EMPTYRECORD Prototype for one plottable timetable record.
   record = struct('source', "", 'kind', "", 'payload', timetable.empty, ...
      'variables', strings(1, 0), 'name', "");
end

function record = makeRecord(source, kind, T)
   %MAKERECORD Build one timetable record with numeric variable inventory.
   variables = numericVariables(T);
   record = struct('source', source, 'kind', kind, 'payload', T, ...
      'variables', variables, 'name', sourceDisplayName(source, kind));
end

function label = sourceDisplayName(source, kind)
   %SOURCEDISPLAYNAME Return a human provenance label for one artifact record.
   source = lower(strtrim(string(source)));
   switch source
      case "promice"
         label = "PROMICE";
      case "esm_snowmip"
         label = "ESM-SnowMIP";
      case "laugh_tests"
         label = "Laugh Tests";
      case "sumup"
         label = "SUMup";
      case "retmip"
         label = "RetMIP";
      case "imau"
         label = "IMAU";
      case "research_site"
         label = "Research site";
      case "mar3.11"
         label = "MAR 3.11";
      case "merra2"
         label = "MERRA-2";
      case "racmo2.3p3"
         label = "RACMO 2.3p3";
      case "modis"
         label = "MODIS (GEUS)";
      otherwise
         label = string(source);
   end
    if isObservationRole(source, kind)
       % Native evaluation userdata remains observational even when it is
       % loaded from the staged userdata tree rather than a target bundle.
       label = label + " observations";
   elseif kind == "met"
      label = label + " forcing";
   end
end

function tf = isObservationRole(source, kind)
   %ISOBSERVATIONROLE True for target data and native evaluation userdata.

   native_sources = ["promice", "retmip", "imau", "esm_snowmip", ...
      "sumup", "research_site"];
   tf = lower(string(kind)) == "observations" ...
      || (lower(string(kind)) == "userdata" ...
      && ismember(lower(string(source)), native_sources));
end

function T = tableToTimeRecord(T)
   %TABLETOTIMERECORD Convert supported time-column tables to timetables.
   names = string(T.Properties.VariableNames);
   time_name = names(find(ismember(lower(names), ["time", "datetime", ...
      "date"]), 1));
   if ~isempty(time_name)
      T = table2timetable(T, 'RowTimes', char(time_name));
      return
   end

   if all(ismember(["start_date", "end_date"], names))
      T = intervalTableToTimetable(T, names);
      return
   end

   T = [];
end

function TT = intervalTableToTimetable(T, names)
   %INTERVALTABLETOTIMETABLE Represent accumulated observations by their spans.
   start_time = T.start_date;
   end_time = T.end_date;
   if ~isdatetime(start_time)
      start_time = datetime(start_time);
   end
   if ~isdatetime(end_time)
      end_time = datetime(end_time);
   end
   % Verification artifacts use UTC. Assign UTC to naive calendar labels and
   % convert already zoned datetimes so interval lines can share model axes.
   start_time.TimeZone = 'UTC';
   end_time.TimeZone = 'UTC';

   value_names = names(arrayfun(@(name) isnumeric(T.(char(name))) ...
      && isvector(T.(char(name))) && ~isStructuralVariable(name), names));
   if isempty(value_names)
      TT = timetable.empty;
      return
   end

   % A third duplicate endpoint with a missing value breaks the line before
   % the next independent observation interval.
   time = reshape([start_time(:) end_time(:) end_time(:)].', [], 1);
   TT = timetable('RowTimes', time);
   for name = reshape(value_names, 1, [])
      values = T.(char(name));
      TT.(char(name)) = reshape( ...
         [values(:) values(:) nan(size(values(:)))].', [], 1);
   end
   TT.Properties.UserData = struct( ...
      'interval_observation', true, ...
      'interval_variables', value_names(:));
   TT = copyTableVariableMetadata(TT, T, value_names);
end

function TT = copyTableVariableMetadata(TT, T, value_names)
   %COPYTABLEVARIABLEMETADATA Preserve units/descriptions on interval records.
   old_names = string(T.Properties.VariableNames);
   new_names = string(TT.Properties.VariableNames);
   for name = reshape(value_names, 1, [])
      old_idx = find(old_names == name, 1);
      new_idx = find(new_names == name, 1);
      if isempty(old_idx) || isempty(new_idx)
         continue
      end
      if numel(T.Properties.VariableUnits) >= old_idx
         TT.Properties.VariableUnits{new_idx} = ...
            T.Properties.VariableUnits{old_idx};
      end
      if numel(T.Properties.VariableDescriptions) >= old_idx
         TT.Properties.VariableDescriptions{new_idx} = ...
            T.Properties.VariableDescriptions{old_idx};
      end
   end
end

function variables = numericVariables(T)
   %NUMERICVARIABLES Return numeric vector variables from one timetable.
   names = string(T.Properties.VariableNames);
   keep = false(1, numel(names));
   for k = 1:numel(names)
      value = T.(char(names(k)));
      keep(k) = isnumeric(value) && isvector(value) ...
         && ~isStructuralVariable(names(k));
   end
   variables = names(keep);
end

function tf = isStructuralVariable(name)
   %ISSTRUCTURALVARIABLE True for coordinates, keys, axes, and support metadata.
   structural = ["time", "datetime", "date", "latitude", "longitude", ...
      "elevation", "elev", "depth", "depth_m", "start_depth", ...
      "stop_depth", "midpoint", "duration", "start_year", "end_year", ...
      "timestamp", "name_key", "reference_key", "method_key", ...
      "measurement_id", "error"];
   tf = ismember(lower(name), structural);
end

function profiles = profileRecords(target)
   %PROFILERECORDS Collect depth-profile tables from observations.
   payload = target.payload;
   profile_cells = cell(1, 4);
   n_cells = 0;
   if isfield(payload, 'data')
      n_cells = n_cells + 1;
      profile_cells{n_cells} = profilesFromValue(payload.data, ...
         "observations");
   end
   for name = ["density", "subsurface_temperature", "smb"]
      if isfield(payload, char(name))
         n_cells = n_cells + 1;
         profile_cells{n_cells} = profilesFromValue( ...
            payload.(char(name)), name);
      end
   end
   profiles = concatProfiles(profile_cells(1:n_cells));
end

function profiles = profilesFromValue(value, source)
   %PROFILESFROMVALUE Recursively collect tables that carry a depth axis.
   profiles = emptyProfile();
   profiles(1) = [];
   if istable(value) || istimetable(value)
      names = string(value.Properties.VariableNames);
      depth_vars = ["depth_m", "depth", "midpoint", "start_depth"];
      if any(ismember(depth_vars, names))
         vars = numericProfileVariables(value, depth_vars);
         if ~isempty(vars)
            profiles = struct('source', source, 'payload', value, ...
               'value_variables', vars);
         end
      end
   elseif isstruct(value)
      fields = fieldnames(value);
      profile_cells = cell(1, numel(fields));
      for k = 1:numel(fields)
         profile_cells{k} = profilesFromValue(value.(fields{k}), ...
            source + ":" + string(fields{k}));
      end
      profiles = concatProfiles(profile_cells);
   end
end

function profiles = concatProfiles(profile_cells)
   %CONCATPROFILES Concatenate optional profile arrays from recursive readers.
   profiles = emptyProfile();
   profiles(1) = [];
   if isempty(profile_cells)
      return
   end
   keep = ~cellfun(@isempty, profile_cells);
   profile_cells = profile_cells(keep);
   if ~isempty(profile_cells)
      profiles = [profile_cells{:}];
   end
end

function profile = emptyProfile()
   %EMPTYPROFILE Prototype for one depth-profile record.
   profile = struct('source', "", 'payload', table(), ...
      'value_variables', strings(1, 0));
end

function vars = numericProfileVariables(T, depth_vars)
   %NUMERICPROFILEVARIABLES Return finite numeric profile value variables.
   names = string(T.Properties.VariableNames);
   depth_name = depth_vars(find(ismember(depth_vars, names), 1));
   depth = T.(char(depth_name));
   vars = strings(1, numel(names));
   n_vars = 0;
   for k = 1:numel(names)
      if ismember(names(k), depth_vars) || isStructuralVariable(names(k))
         continue
      end
      value = T.(char(names(k)));
      if isnumeric(value) && isvector(value) ...
            && any(isfinite(value) & isfinite(depth))
         n_vars = n_vars + 1;
         vars(n_vars) = names(k);
      end
   end
   vars = vars(1:n_vars);
end

function [payloads, names, observations, albedo_sources] = ...
      recordsForVariable(records, varname)
   %RECORDSFORVARIABLE Return payloads, labels, and source roles.
   keep = arrayfun(@(r) ismember(varname, r.variables), records);
   selected = preferObservationRecords(records(keep));
   if varname == "albedo"
      modis = canonicalModisRecord(records);
      if ~isempty(modis)
         selected = [selected, modis];
      end
   end
   payloads = cell(1, numel(selected));
   names = strings(1, numel(selected));
   observations = false(1, numel(selected));
   albedo_sources = repmat("radiometer", 1, numel(selected));
   state_sources = icemodel.verification.namelists.rcmProductIds( ...
      ["mar", "merra"]);
   ratio_sources = icemodel.verification.namelists.rcmProductIds("racmo");
   for k = 1:numel(selected)
      payloads{k} = selected(k).payload;
      names(k) = selected(k).name;
      observations(k) = isObservationRole( ...
         selected(k).source, selected(k).kind);
      if ismember(selected(k).source, state_sources)
         % MAR and MERRA-2 expose albedo as native model state variables.
         albedo_sources(k) = "model_state";
      elseif ismember(selected(k).source, ratio_sources)
         % RACMO albedo is derived from its staged net/incoming shortwave ratio.
         albedo_sources(k) = "model_ratio";
      end
   end
end

function selected = preferObservationRecords(selected)
   %PREFEROBSERVATIONRECORDS Drop same-source userdata when the target has VARNAME.
   if isempty(selected)
      return
   end
   sources = unique([selected.source], 'stable');
   keep = true(size(selected));
   for source = sources
      same_source = [selected.source] == source;
      if any(same_source & [selected.kind] == "observations")
         keep(same_source & [selected.kind] ~= "observations") = false;
      end
   end
   selected = selected(keep);
end

function record = canonicalModisRecord(records)
   %CANONICALMODISRECORD Return one longest-coverage GEUS MODIS albedo series.
   record = emptyRecord();
   record(1) = [];
   keep = arrayfun(@(r) ismember("modis", r.variables), records);
   candidates = records(keep);
   if isempty(candidates)
      return
   end

   coverage = zeros(1, numel(candidates));
   for k = 1:numel(candidates)
      coverage(k) = nnz(isfinite(candidates(k).payload.modis));
   end
   [~, best] = max(coverage);
   payload = candidates(best).payload(:, "modis");
   payload = renamevars(payload, "modis", "albedo");
   record = makeRecord("modis", "external", payload);
end

function [tables, names] = profilesForVariable(profiles, varname, dataset_family)
   %PROFILESFORVARIABLE Return human-labelled profile series carrying VARNAME.
   keep = arrayfun(@(p) ismember(varname, p.value_variables), profiles);
   selected = profiles(keep);
   tables = cell(1, 0);
   names = strings(1, 0);
   base_name = sourceDisplayName(dataset_family, "observations");
   for k = 1:numel(selected)
      [next_tables, next_names] = splitNamedProfiles( ...
         selected(k).payload, base_name);
      tables = [tables, next_tables]; %#ok<AGROW>
      names = [names, next_names]; %#ok<AGROW>
   end
end

function [tables, names] = splitNamedProfiles(T, base_name)
   %SPLITNAMEDPROFILES Color a small number of explicit source series separately.

   tables = {T};
   names = base_name;
   if ~ismember("name", string(T.Properties.VariableNames))
      return
   end

   labels = string(T.name);
   labelled = ~ismissing(labels) & strlength(strtrim(labels)) > 0;
   series = unique(labels(labelled), 'stable');
   if numel(series) < 2 || numel(series) > 24
      % Large or absent identity sets would produce an unreadable legend; in
      % that case retain the source-faithful aggregate profile collection.
      return
   end

   tables = cell(1, numel(series) + any(~labelled));
   names = strings(1, numel(tables));
   for k = 1:numel(series)
      tables{k} = T(labels == series(k), :);
      names(k) = base_name + " - " + series(k);
   end
   if any(~labelled)
      tables{end} = T(~labelled, :);
      names(end) = base_name + " - unlabeled";
   end
end

function vars = allNumericVariables(records, profiles)
   %ALLNUMERICVARIABLES Union all numeric timetable and profile variables.
   vars = strings(1, 0);
   if ~isempty(records)
      vars = unique([vars, records.variables], 'stable');
   end
   if ~isempty(profiles)
      vars = unique([vars, profiles.value_variables], 'stable');
   end
end

function noDataPanel(ax, varname, reason)
   %NODATAPANEL Render a visible no-data tile.
   title(ax, strrep(varname, '_', '\_'), 'FontWeight', 'bold')
   text(ax, 0.5, 0.5, "no plottable data" + newline + "(" + reason + ")", ...
      'Units', 'normalized', 'HorizontalAlignment', 'center', ...
      'VerticalAlignment', 'middle', 'Color', [0.35 0.35 0.35])
   ax.XTick = [];
   ax.YTick = [];
   ax.Box = 'on';
end

function shareFiniteTimeExtent(panel_axes)
   %SHAREFINITETIMEEXTENT Apply one data-supported time range to every panel.
   first_time = NaT;
   last_time = NaT;
   time_axes = gobjects(numel(panel_axes), 1);
   n_time_axes = 0;

   % Use finite plotted samples rather than timetable endpoints so missing edge
   % values do not extend a figure beyond the support visible in any panel.
   for ax = reshape(panel_axes, 1, [])
      lines = findall(ax, 'Type', 'Line');
      supports_time = false;
      for line = reshape(lines, 1, [])
         x = line.XData;
         y = line.YData;
         if ~isdatetime(x) || ~isnumeric(y) || numel(x) ~= numel(y)
            continue
         end
         finite = ~isnat(x) & isfinite(y);
         if ~any(finite)
            continue
         end
         supports_time = true;
         line_first = min(x(finite));
         line_last = max(x(finite));
         if isnat(first_time) || line_first < first_time
            first_time = line_first;
         end
         if isnat(last_time) || line_last > last_time
            last_time = line_last;
         end
      end
      if supports_time
         n_time_axes = n_time_axes + 1;
         time_axes(n_time_axes) = ax;
      end
   end

   if n_time_axes == 0
      return
   end
   time_axes = time_axes(1:n_time_axes);
   if first_time == last_time
      % A one-sample figure still needs increasing limits for MATLAB's ruler.
      first_time = first_time - days(0.5);
      last_time = last_time + days(0.5);
   end

   % Link after finding the union so later interactive zooming also remains
   % synchronized across the complete tiled diagnostic figure.
   linkaxes(time_axes, 'x')
   xlim(time_axes(1), [first_time, last_time])
end

function lgd = styleLegend(ax, reserve_space)
   %STYLELEGEND Keep compact legends readable under light and dark themes.
   if nargin < 2
      reserve_space = true;
   end
   lgd = gobjects(0, 1);
   handles = findobj(ax, '-property', 'DisplayName');
   % FINDOBJ can return heterogeneous graphics subclasses. GET performs the
   % shared-property dispatch safely where comma-separated dot expansion does
   % not (for example, axes containing internal scribe peers).
   display_names = get(handles, 'DisplayName');
   if ~iscell(display_names)
      display_names = {display_names};
   end
   handles = handles(cellfun( ...
      @(name) any(strlength(string(name)) > 0, 'all'), display_names));
   if isempty(handles)
      return
   end

   % Graphics searches return reverse creation order; restore plot order so
   % source and paired-variable labels remain stable across every subset.
   handles = flipud(handles(:));
   n_columns = min(4, numel(handles));
   lgd = legend(ax, handles, 'Location', 'best', 'Interpreter', 'none', ...
      'NumColumns', n_columns, ...
      'FontSize', 7, 'Color', 'w', 'TextColor', 'k', ...
      'EdgeColor', [0.75 0.75 0.75]);
   if reserve_space
      % A small symmetric seed margin gives BEST enough room to choose a stable
      % edge before the final one-sided in-axes clearance pass. Specialized
      % legends moved outside immediately after creation skip this seed.
      padYLimits(ax, 0.03, 0.03)
   end
end

function finalizeLegendClearance(fig)
   %FINALIZELEGENDCLEARANCE Reserve space after tiled geometry has settled.

   % Best-position legends can move while later tiles and shared time limits
   % are added. Resolve their owning axes by final pixel geometry immediately
   % before export so one shared pass protects every stacked figure family.
   drawnow
   axes_handles = findall(fig, 'Type', 'Axes');
   legends = findall(fig, 'Type', 'Legend');
   for lgd = reshape(legends, 1, [])
      if contains(string(lgd.Location), "outside")
         continue
      end
      legend_pixels = getpixelposition(lgd, true);
      legend_center = legend_pixels(1:2) + legend_pixels(3:4) / 2;
      owner = gobjects(0, 1);
      for ax = reshape(axes_handles, 1, [])
         axes_pixels = getpixelposition(ax, true);
         if legend_center(1) >= axes_pixels(1) ...
               && legend_center(1) <= axes_pixels(1) + axes_pixels(3) ...
               && legend_center(2) >= axes_pixels(2) ...
               && legend_center(2) <= axes_pixels(2) + axes_pixels(4)
            owner(end + 1, 1) = ax; %#ok<AGROW>
         end
      end
      if isscalar(owner)
         % Freeze MATLAB's resolved best position before changing limits;
         % otherwise BEST can move the legend after clearance is computed.
         resolved_position = lgd.Position;
         lgd.Location = 'none';
         lgd.Position = resolved_position;
         reserveLegendClearance(owner, lgd, 0.03)
      end
   end
end

function reserveLegendClearance(ax, lgd, base_fraction)
   %RESERVELEGENDCLEARANCE Protect extrema beneath an in-axes legend.

   % Layout geometry is finalized lazily. LIMITRATE avoids executing callbacks
   % while still resolving tiled axes and legend pixel positions for export.
   drawnow limitrate nocallbacks
   axes_pixels = getpixelposition(ax, true);
   legend_pixels = getpixelposition(lgd, true);
   if any(~isfinite([axes_pixels, legend_pixels])) || axes_pixels(4) <= 0
      padYLimits(ax, base_fraction, base_fraction)
      return
   end

   % Measure the legend against the axes plot box. The axes pixel rectangle
   % includes a small inset around the plotted box, so retain a four-percent
   % visual gap between line crowns and the legend border. Middle legends keep
   % the small symmetric fallback because changing either limit cannot clear a
   % genuinely central overlay.
   relative_bottom = (legend_pixels(2) - axes_pixels(2)) / axes_pixels(4);
   relative_top = (legend_pixels(2) + legend_pixels(4) - axes_pixels(2)) ...
      / axes_pixels(4);
   legend_center = (relative_bottom + relative_top) / 2;
   geometry_gap = 0.04;
   visual_lower = base_fraction;
   visual_upper = base_fraction;
   if legend_center >= 2 / 3
      visual_lower = 0;
      visual_upper = max(base_fraction, 1 - relative_bottom + geometry_gap);
   elseif legend_center <= 1 / 3
      visual_lower = max(base_fraction, relative_top + geometry_gap);
      visual_upper = 0;
   end

   % Compact legends should never consume more than one quarter of the data
   % range. The cap protects unusual multirow legends while the four-column
   % policy keeps normal stacked-panel legends well below it.
   visual_lower = min(0.24, max(0, visual_lower));
   visual_upper = min(0.24, max(0, visual_upper));
   if string(ax.YDir) == "reverse"
      padYLimits(ax, visual_upper, visual_lower)
   else
      padYLimits(ax, visual_lower, visual_upper)
   end
end

function padYLimits(ax, lower_fraction, upper_fraction)
   %PADYLIMITS Reserve lower and upper fractions of the final y-axis range.
   for ruler = reshape(ax.YAxis, 1, [])
      limits = ruler.Limits;
      span = diff(limits);
      if isfinite(span) && span > 0
         % Work from MATLAB's resolved limits so constant and sparse series keep
         % the native safe range selected by the graphics engine. Solving for
         % the final span makes each requested fraction a true visual reserve.
         final_span = span / (1 - lower_fraction - upper_fraction);
         ruler.Limits = [limits(1) - lower_fraction * final_span, ...
            limits(2) + upper_fraction * final_span];
      end
   end
end

function figfile = exportFigure(fig, case_dir, group_name, save_figs)
   %EXPORTFIGURE Write a figure to disk when requested.
   finalizeLegendClearance(fig)
   figfile = "";
   if ~save_figs
      return
   end
   figfile = fullfile(case_dir, char(group_name + ".png"));
   exportgraphics(fig, figfile, 'Resolution', 150);
end

function setFigureSize(f, width, height)
   %SETFIGURESIZE Use stable pixel sizes for exported verification figures.
   f.Units = 'pixels';
   f.Position(3:4) = [width height];
end

function name = figureName(c, suffix)
   %FIGURENAME Stable figure window name.
   name = sprintf('%s/%s %s', c.dataset_family, c.case_id, suffix);
end

function state = visibleState(visible)
   %VISIBLESTATE MATLAB figure Visible property value.
   state = 'off';
   if visible
      state = 'on';
   end
end

function row = summaryRow(c, group_name, figfile, plotted, unplotted, ...
      target, met, data)
   %SUMMARYROW Build one output table row.
   period = manifestPeriod(c);
   row = table(string(c.dataset_family), string(c.case_id), ...
      string(group_name), string(target.format), ...
      string(strjoin(target.variables, ",")), ...
      string(strjoin(met.sources, ",")), ...
      string(strjoin(data.sources, ",")), ...
      string(strjoin(plotted, ",")), string(strjoin(unplotted, ",")), ...
      string(period.start), string(period.end), string(figfile), ...
      'VariableNames', {'dataset_family', 'case_id', 'figure_group', ...
      'target_format', 'target_variables', 'met_sources', ...
      'userdata_sources', 'plotted_variables', 'unplotted_variables', ...
      'period_start', 'period_end', 'figure_file'});
end

function row = addUnplotted(row, missing)
   %ADDUNPLOTTED Store case-level unplotted variables on the final row.
   row.unplotted_variables = string(strjoin(missing, ","));
end

function period = manifestPeriod(c)
   %MANIFESTPERIOD Return a stable start/end struct from a manifest row.
   period = icemodel.verification.helpers.fieldOr(c, 'period', ...
      struct('start', '', 'end', ''));
   if ~isfield(period, 'start')
      period.start = '';
   end
   if ~isfield(period, 'end')
      period.end = '';
   end
end

function summary = emptySummary()
   %EMPTYSUMMARY Stable empty summary table for no manifest/case matches.
   summary = table(strings(0, 1), strings(0, 1), strings(0, 1), ...
      strings(0, 1), strings(0, 1), strings(0, 1), strings(0, 1), ...
      strings(0, 1), strings(0, 1), strings(0, 1), strings(0, 1), ...
      strings(0, 1), ...
      'VariableNames', {'dataset_family', 'case_id', 'figure_group', ...
      'target_format', 'target_variables', 'met_sources', ...
      'userdata_sources', 'plotted_variables', 'unplotted_variables', ...
      'period_start', 'period_end', 'figure_file'});
end
