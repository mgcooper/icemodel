function mustBeUservarName(arg)
%MUSTBEUSERVARNAME Validate that input is a valid userdata variable name.

   if isempty(arg)
      return
   end

   valid = icemodel.namelists.uservars();
   if ~ismember(string(arg), valid)
      eid = ['icemodel:validators:' mfilename];
      msg = sprintf('Value must be one of:\n%s', strjoin(cellstr(valid), ', '));
      throwAsCaller(MException(eid, msg));
   end
end
