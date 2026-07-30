function out = profile(data, value_name, kwargs)
   %PROFILE Plot one or more depth-profile tables.
   %
   %  out = icemodel.plot.profile({density1, density2}, "density_kgm3")
   %
   % The helper keeps verification profile plots consistent: value on x, depth
   % on y, reversed y-axis, and source names in the legend.

   arguments
      data
      value_name (1, 1) string
      kwargs.depth_name (1, 1) string = ""
      kwargs.names string = strings(0, 1)
      kwargs.axes matlab.graphics.axis.Axes = gca
      kwargs.startdate = ""
      kwargs.enddate = ""
   end

   tables = normalizeTables(data);
   names = normalizeNames(kwargs.names, numel(tables));
   colors = lines(max(numel(tables), 1));
   handles = gobjects(numel(tables), 1);
   plotted = false(numel(tables), 1);
   x_unit = "";
   y_unit = "";

   was_hold = ishold(kwargs.axes);
   hold(kwargs.axes, 'on')
   for k = 1:numel(tables)
      T = icemodel.plot.filterDateWindow( ...
         tables{k}, kwargs.startdate, kwargs.enddate);
      if isempty(T)
         continue
      end
      depth_name = resolveDepthName(T, kwargs.depth_name);
      if depth_name == "" || ~ismember(value_name, ...
            string(T.Properties.VariableNames))
         continue
      end
      x = T.(char(value_name));
      y = T.(char(depth_name));
      if ~isnumeric(x) || ~isnumeric(y)
         continue
      end
      if x_unit == ""
         x_unit = icemodel.plot.variableUnit(T, value_name);
      end
      if y_unit == ""
         y_unit = icemodel.plot.variableUnit(T, depth_name);
      end
      handles(k) = plot(kwargs.axes, x, y, '.', 'Color', colors(k, :), ...
         'DisplayName', char(names(k)));
      plotted(k) = true;
   end
   if ~was_hold
      hold(kwargs.axes, 'off')
   end

   set(kwargs.axes, 'YDir', 'reverse')
   xlabel(kwargs.axes, labelText(value_name, x_unit), 'Interpreter', 'none')
   ylabel(kwargs.axes, labelText("depth", y_unit), 'Interpreter', 'none')
   grid(kwargs.axes, 'on')
   if any(plotted)
      legend(kwargs.axes, 'Location', 'best', 'Interpreter', 'none')
   end

   out = struct('axes', kwargs.axes, 'lines', handles(plotted), ...
      'plotted', plotted);
end

function tables = normalizeTables(data)
   %NORMALIZETABLES Convert a table or cell array to a row cell array.

   if istable(data) || istimetable(data)
      if istimetable(data)
         data = icemodel.plot.canonicalTimeDimension(data);
      end
      tables = {data};
      return
   end
   if iscell(data)
      tables = reshape(data, 1, []);
      for k = 1:numel(tables)
         if istimetable(tables{k})
            tables{k} = icemodel.plot.canonicalTimeDimension(tables{k});
         end
      end
      return
   end
   error('icemodel:plot:profile:badInput', ...
      'data must be a table/timetable or cell array of tables')
end

function names = normalizeNames(names, n_tables)
   %NORMALIZENAMES Return one display name per profile table.

   if isempty(names)
      names = "profile " + string(1:n_tables);
      return
   end
   names = reshape(names, 1, []);
   if isscalar(names)
      names = repmat(names, 1, n_tables);
   elseif numel(names) ~= n_tables
      error('icemodel:plot:profile:badNames', ...
         'names must be scalar or match the number of tables')
   end
end

function depth_name = resolveDepthName(T, requested)
   %RESOLVEDEPTHNAME Find the depth axis in a profile table.

   if requested ~= ""
      depth_name = requested;
      return
   end
   names = string(T.Properties.VariableNames);
   hits = ["depth_m", "depth", "midpoint", "start_depth"];
   depth_name = hits(find(ismember(hits, names), 1));
   if isempty(depth_name)
      depth_name = "";
   end
end

function label = labelText(varname, unit)
   %LABELTEXT Compose a compact variable + unit label.

   label = varname;
   if unit ~= ""
      label = varname + " [" + unit + "]";
   end
end
