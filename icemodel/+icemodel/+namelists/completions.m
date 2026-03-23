function list = completions()
%COMPLETIONS Return the supported icemodel.completions selector names.
%
%  list = icemodel.namelists.completions()

   tmp = dir(fullfile(icemodel.internal.fullpath, 'icemodel', '+icemodel', ...
      '+namelists', '*.m'));
   list = string(sort(strrep({tmp.name}.', '.m', '')));
end
