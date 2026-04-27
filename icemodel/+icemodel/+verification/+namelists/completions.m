function list = completions()
   %COMPLETIONS Return the supported verification namelist selector names.
   %
   %  list = icemodel.verification.namelists.completions()
   %
   % Outputs
   %  list   String column of namelist function names.
   %
   % Role
   %  Convenience catalog for interactive discovery of verification namelists.

   % Discover namelist files from disk so this catalog stays current when new
   % selector functions are added.
   files = dir(fullfile(icemodel.internal.fullpath('icemodel'), ...
      '+icemodel', '+verification', '+namelists', '*.m'));
   names = string(erase({files.name}, '.m'))';
   list = sort(names); %#ok<TRSRT>
end
