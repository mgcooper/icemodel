function value = fieldOr(s, fieldname, default_value)
   %FIELDOR Return a struct field or default value.
   %
   %  value = icemodel.verification.helpers.fieldOr(s, fieldname, default_value)
   %
   % Manifest readers use this to tolerate additive schema fields without
   % scattering isfield/default logic across plotting and setup code.

   if isstruct(s) && isfield(s, fieldname)
      value = s.(fieldname);
   else
      value = default_value;
   end
end
