function setup()
   thispath = fileparts(mfilename('fullpath'));
   addpath(genpath(thispath))
   icemodel.config('casename', 'demo');
end
