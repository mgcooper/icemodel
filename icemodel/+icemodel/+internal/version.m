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
   % citation file. The override and reset branches below update the cache.
   persistent current
   is_reset = nargin == 1 && strcmp('reset', new);
   if nargin == 1 && ~is_reset && isrow(new) && ischar(new)
      % Apply an override before any cache check, so an override is not lost
      % on the first call after the cache is cleared.
      current = new;
   elseif isempty(current) || is_reset
      % Discard an override before reading so a failed reset cannot leave
      % stale process state in place of a missing or malformed version source.
      current = [];
      % Resolve from this installed package rather than the fullpath helper,
      % so version lookup stays within the documented MATLAB R2017a floor.
      internal_dir = fileparts(mfilename('fullpath'));
      package_dir = fileparts(internal_dir);
      toolbox_dir = fileparts(package_dir);
      project_dir = fileparts(toolbox_dir);
      current = icemodel.internal.readCffVersion( ...
         fullfile(project_dir, 'CITATION.cff'));
   end
   str = current;
end
