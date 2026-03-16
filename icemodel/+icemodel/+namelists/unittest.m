function list = unittest()
%UNITTEST Return the supported unit test file names.
%
%  list = icemodel.namelists.unittest()

   tmp = dir(fullfile(icemodel.internal.fullpath, 'test', 'unit', '*.m'));
   list = string(sort(erase({tmp.name}.', '.m')));
end
