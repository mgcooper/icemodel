function mustBePromiceImportSourceSelection(arg)
   %MUSTBEPROMICEIMPORTSOURCESELECTION Validate PROMICE import source selectors.
   %
   % PROMICE imports accept the native PROMICE source plus the shared verification
   % RCM sources.
   values = reshape(string(arg), 1, []);
   values = values(strlength(strtrim(values)) > 0);
   if isempty(values)
      return
   end

   valid = ["promice", icemodel.verification.namelists.rcmsources()];
   bad = setdiff(values, valid, 'stable');
   if ~isempty(bad)
      eid = ['icemodel:validators:' mfilename];
      msg = sprintf('Value must be one of:\n%s', strjoin(cellstr(valid), ', '));
      throwAsCaller(MException(eid, msg));
   end
end
