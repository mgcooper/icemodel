function setup()
   thispath = fileparts(mfilename('fullpath'));
   addpath(genpath(thispath))

   setenv('ICEMODEL_INPUT_PATH', fullfile(thispath, 'demo', 'input'))
   setenv('ICEMODEL_OUTPUT_PATH', fullfile(thispath, 'demo', 'output'))
end
