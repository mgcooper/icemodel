function tbl = renameTableVars(tbl, oldvars, newvars)
   %RENAMETABLEVARS Rename variables in a table.
   %
   %   TBL = renameTableVars(TBL, OLDVARS, NEWVARS) renames the variables in the
   %   table TBL. OLDVARS is a cell array or string array of the current variable
   %   names, and NEWVARS is a cell array or string array of the new variable
   %   names. OLDVARS and NEWVARS must have the same length.
   %
   %   Note that this function is deliberately designed for pre-R2020a
   %   compatibility. For MATLAB versions >=R2020a, use renamevars.
   %
   %   Example:
   %       T = table([1; 2], [3; 4], 'VariableNames', {'Var1', 'Var2'});
   %       T = renameTableVars(T, {'Var1', 'Var2'}, {'NewVar1', 'NewVar2'});
   %
   %   Input Arguments:
   %       TBL     - Input table
   %       OLDVARS - Cell array or string array of old variable names
   %       NEWVARS - Cell array or string array of new variable names
   %
   %   Output Arguments:
   %       TBL     - Output table with renamed variables
   %
   %   See also: renamevars

   % Ensure oldvars and newvars are cell arrays of strings
   if ischar(oldvars)
      oldvars = {oldvars};
   end
   if ischar(newvars)
      newvars = {newvars};
   end

   % Validate attributes of input arguments
   validateattributes(oldvars, ...
      {'cell', 'string'}, {'nonempty'}, mfilename, 'OLDVARS', 2);
   validateattributes(newvars, ...
      {'cell', 'string'}, {'nonempty'}, mfilename, 'NEWVARS', 3);

   % Convert to string arrays for consistency
   oldvars = string(oldvars);
   newvars = string(newvars);

   % Ensure oldvars and newvars have the same length
   assert(numel(oldvars) == numel(newvars), ...
      'OLDVARS and NEWVARS must have the same length.');

   % Get the current variable names of the table
   varnames = string(tbl.Properties.VariableNames);

   % Loop through oldvars and replace with newvars
   for n = 1:numel(oldvars)
      if any(strcmp(oldvars(n), varnames))
         tbl.Properties.VariableNames(oldvars(n)) = newvars(n);
      else
         warning('Variable name %s not found in the table.', oldvars(n));
      end
   end
end
