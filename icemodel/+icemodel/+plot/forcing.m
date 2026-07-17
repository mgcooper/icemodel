function out = forcing(inputs, kwargs)
   %FORCING Plot one or more icemodel met-file forcings.
   %
   %  out = icemodel.plot.forcing("met_kanm_mar_2012_1hr.mat")
   %  out = icemodel.plot.forcing({met1, met2}, names=["mar","merra"])
   %
   % Inputs can be paths to MAT files containing `met`, timetables, or cell/
   % string arrays of either. Energy fluxes are retimed to daily means by
   % default and plotted on the left axis; near-surface air temperature is
   % plotted in degC on the right axis for context.

   arguments
      inputs
      kwargs.names string = strings(0, 1)
      kwargs.frequency (1, 1) string ...
         {mustBeMember(kwargs.frequency, ["hourly", "daily"])} = "daily"
      kwargs.axes matlab.graphics.axis.Axes = gca
      kwargs.title (1, 1) string = "staged met files"
      kwargs.startdate = ""
      kwargs.enddate = ""
   end

   [tables, names] = loadInputs(inputs, kwargs.names);
   tables = filterInputs(tables, kwargs.startdate, kwargs.enddate);
   energy_vars = ["swd", "swu", "swn", "lwd", "lwu", "lwn", ...
      "netr", "shf", "lhf", "thf"];
   colors = icemodel.plot.sourceColor(names);
   energy_handles = cell(numel(tables) * numel(energy_vars), 1);
   tair_handles = cell(numel(tables), 1);
   n_energy = 0;
   n_tair = 0;

   title(kwargs.axes, kwargs.title, 'Interpreter', 'none')
   yyaxis(kwargs.axes, 'left')
   hold(kwargs.axes, 'on')
   for k = 1:numel(tables)
      T = retime(tables{k}, char(kwargs.frequency), 'mean');
      present = energy_vars(ismember(energy_vars, ...
         string(T.Properties.VariableNames)));
      for n = 1:numel(present)
         h = plot(kwargs.axes, T.Time, T.(char(present(n))), ...
            'LineWidth', 1.0, 'Color', colors(k, :), ...
            'DisplayName', char(names(k) + ":" + present(n)));
         n_energy = n_energy + 1;
         energy_handles{n_energy} = h;
      end
   end
   ylabel(kwargs.axes, 'energy flux [W m^{-2}]')

   yyaxis(kwargs.axes, 'right')
   hold(kwargs.axes, 'on')
   for k = 1:numel(tables)
      T = retime(tables{k}, char(kwargs.frequency), 'mean');
      if ~ismember("tair", string(T.Properties.VariableNames))
         continue
      end
      tair = T.tair;
      if median(tair, 'omitnan') > 100
         tair = tair - 273.15;
      end
      h = plot(kwargs.axes, T.Time, tair, '--', 'LineWidth', 1.2, ...
         'Color', colors(k, :), 'DisplayName', char(names(k) + ":tair"));
      n_tair = n_tair + 1;
      tair_handles{n_tair} = h;
   end
   ylabel(kwargs.axes, 'tair [degC]')

   grid(kwargs.axes, 'on')
   energy_handles = handleVector(energy_handles(1:n_energy));
   tair_handles = handleVector(tair_handles(1:n_tair));
   all_handles = [energy_handles; tair_handles];
   if ~isempty(all_handles)
      lgd = legend(kwargs.axes, all_handles, 'Location', 'best', ...
         'Interpreter', 'none', 'FontSize', 7, 'Color', 'w', ...
         'TextColor', 'k', 'EdgeColor', [0.75 0.75 0.75]);
      lgd.NumColumns = min(3, numel(lgd.String));
   end

   out = struct('axes', kwargs.axes, 'energy', energy_handles, ...
      'tair', tair_handles, 'names', names);
end

function tables = filterInputs(tables, startdate, enddate)
   %FILTERINPUTS Apply the same optional date window to every met timetable.

   for k = 1:numel(tables)
      tables{k} = icemodel.plot.filterDateWindow(tables{k}, startdate, enddate);
   end
end

function [tables, names] = loadInputs(inputs, names)
   %LOADINPUTS Normalize input paths/timetables to tables plus display names.

   cells = inputCells(inputs);
   tables = cell(1, numel(cells));
   default_names = strings(1, numel(cells));
   for k = 1:numel(cells)
      [tables{k}, default_names(k)] = loadOne(cells{k});
   end
   if isempty(names)
      names = default_names;
   else
      names = reshape(names, 1, []);
      if isscalar(names)
         names = repmat(names, 1, numel(cells));
      elseif numel(names) ~= numel(cells)
         error('icemodel:plot:forcing:badNames', ...
            'names must be scalar or match the number of inputs')
      end
   end
end

function cells = inputCells(inputs)
   %INPUTCELLS Convert supported input containers to a cell array.

   if istimetable(inputs)
      cells = {inputs};
   elseif iscell(inputs)
      cells = reshape(inputs, 1, []);
   elseif isstring(inputs) || ischar(inputs)
      cells = cellstr(string(inputs));
   else
      error('icemodel:plot:forcing:badInput', ...
         'inputs must be met paths, timetables, or cell/string arrays')
   end
end

function handles = handleVector(cells)
   %HANDLEVECTOR Convert a collected handle cell to a graphics handle column.
   if isempty(cells)
      handles = gobjects(0, 1);
   else
      handles = vertcat(cells{:});
   end
end

function [T, name] = loadOne(input)
   %LOADONE Load one met timetable from memory or a MAT file.

   if istimetable(input)
      T = icemodel.plot.canonicalTimeDimension(input);
      name = "met";
      return
   end

   pathname = string(input);
   loaded = load(pathname, 'met');
   if ~isfield(loaded, 'met') || ~istimetable(loaded.met)
      error('icemodel:plot:forcing:noMet', ...
         '%s does not contain a timetable named met', pathname)
   end
   T = icemodel.plot.canonicalTimeDimension(loaded.met);
   [~, stem] = fileparts(pathname);
   name = string(stem);
end
