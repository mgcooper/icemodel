function str = version(new)
   %VERSION Set or get the IceModel version number.
   %
   %  VERS = ICEMODEL.INTERNAL.VERSION()
   %  VERS = ICEMODEL.INTERNAL.VERSION(NEW_VERSION)
   %  VERS = ICEMODEL.INTERNAL.VERSION('reset')
   %
   %  The default version comes from the repository CITATION.cff file. A
   %  NEW_VERSION remains a process-local override until reset.
   %
   % See also: ICEMODEL.INTERNAL.READCFFVERSION,
   % ICEMODEL.INTERNAL.REFERENCE

   % Cache the persisted version so ordinary runtime calls do not reread the
   % citation file, while preserving the established override/reset contract.
   persistent current
   if isempty(current) || nargin == 1 && strcmp('reset', new)
      % Discard an override before reading so a failed reset cannot leave stale
      % process state masking a missing or malformed sole version source.
      current = [];
      % Resolve from this installed package rather than the newer fullpath
      % helper so version lookup retains the documented MATLAB R2017a floor.
      internal_dir = fileparts(mfilename('fullpath'));
      package_dir = fileparts(internal_dir);
      toolbox_dir = fileparts(package_dir);
      project_dir = fileparts(toolbox_dir);
      current = icemodel.internal.readCffVersion( ...
         fullfile(project_dir, 'CITATION.cff'));
   elseif nargin == 1 && isrow(new) && ischar(new)
      current = new;
   end
   str = current;
end
