function list = benchmark()
%BENCHMARK Return the supported formal benchmark file names.
%
%  list = icemodel.namelists.benchmark()

   tmp = dir(fullfile(icemodel.internal.fullpath, 'test', 'benchmarks', ...
      '*.m'));
   list = string(sort(erase({tmp.name}.', '.m')));
end
