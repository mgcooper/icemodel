function attrs = sourceVariableAttributes(info)
   %SOURCEVARIABLEATTRIBUTES Preserve NetCDF variable attributes by variable.
   %
   %  attrs = icemodel.forcing.helpers.sourceVariableAttributes(info)
   %
   % Converts the ncinfo variable attribute blocks to a struct keyed by valid
   % MATLAB variable names so builders can preserve source metadata consistently.

   attrs = struct();
   for k = 1:numel(info.Variables)
      v = info.Variables(k);
      one = struct();
      for a = 1:numel(v.Attributes)
         key = matlab.lang.makeValidName(v.Attributes(a).Name);
         one.(key) = string(v.Attributes(a).Value);
      end
      attrs.(matlab.lang.makeValidName(v.Name)) = one;
   end
end
