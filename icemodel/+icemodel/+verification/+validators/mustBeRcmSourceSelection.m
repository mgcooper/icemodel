function mustBeRcmSourceSelection(arg)
   %MUSTBERCMSOURCESELECTION Validate verification RCM source selectors.
   %
   %  icemodel.verification.validators.mustBeRcmSourceSelection(arg)
   %
   % The canonical source labels and order live in
   % icemodel.verification.namelists.rcmsources.
   values = reshape(string(arg), 1, []);
   if isempty(values)
      return
   end

   valid = icemodel.verification.namelists.rcmsources();
   bad = setdiff(values, valid, 'stable');
   if ~isempty(bad)
      eid = ['icemodel:validators:' mfilename];
      msg = sprintf('Value must be one of:\n%s', strjoin(cellstr(valid), ', '));
      throwAsCaller(MException(eid, msg));
   end
end
