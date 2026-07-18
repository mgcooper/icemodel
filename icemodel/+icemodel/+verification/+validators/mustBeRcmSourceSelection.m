function mustBeRcmSourceSelection(arg, native_sources)
   %MUSTBERCMSOURCESELECTION Validate verification forcing-source selectors.
   %
   %  icemodel.verification.validators.mustBeRcmSourceSelection(arg)
   %  icemodel.verification.validators.mustBeRcmSourceSelection( ...
   %     arg, native_sources)
   %
   % The canonical RCM labels and order live in
   % icemodel.verification.namelists.rcmsources. Importers that also expose a
   % native runtime source pass that narrow family source in native_sources.

   arguments
      arg
      native_sources (1, :) string = strings(1, 0)
   end

   values = reshape(string(arg), 1, []);
   values = values(strlength(strtrim(values)) > 0);
   if isempty(values)
      return
   end

   valid = [reshape(native_sources, 1, []), ...
      icemodel.verification.namelists.rcmsources()];
   bad = setdiff(values, valid, 'stable');
   if ~isempty(bad)
      eid = ['icemodel:validators:' mfilename];
      msg = sprintf('Value must be one of:\n%s', strjoin(cellstr(valid), ', '));
      throwAsCaller(MException(eid, msg));
   end
end
