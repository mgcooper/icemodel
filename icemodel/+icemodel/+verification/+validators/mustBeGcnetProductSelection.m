function mustBeGcnetProductSelection(arg)
   %MUSTBEGCNETPRODUCTSELECTION Validate Vandecrux/GC-Net product selectors.
   %
   %  icemodel.verification.validators.mustBeGcnetProductSelection(arg)
   %
   % The canonical selector order lives in gcnetProductNames. Keeping registry
   % lookup inside this validator also makes it safe in MATLAB arguments blocks.

   values = reshape(string(arg), 1, []);
   if isempty(values)
      return
   end

   % Reject unknown selectors with the same stable caller-facing validator id.
   valid = icemodel.verification.setup.gcnetProductNames();
   bad = setdiff(values, valid, 'stable');
   if ~isempty(bad)
      eid = ['icemodel:validators:' mfilename];
      msg = sprintf('Value must be one of:\n%s', strjoin(cellstr(valid), ', '));
      throwAsCaller(MException(eid, msg));
   end
end
