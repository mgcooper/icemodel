function [input_path, output_path, eval_path] = configureModelPaths(rootdir)
%CONFIGUREMODELPATHS Configure icemodel test/build paths.
%
%  [input_path, output_path, eval_path] = test.helpers.configureModelPaths(rootdir)
%
% Formal tests use one fixed public input tree under demo/data/input. This keeps the
% suite independent of the sibling runoff project once the required test files
% are shipped with icemodel.

   % Ensure the repo and test helpers are visible before configuring paths.
   addpath(genpath(rootdir))

   data_path = fullfile(rootdir, 'demo', 'data');
   input_path = fullfile(data_path, 'input');
   output_path = fullfile(data_path, 'output');
   eval_path = fullfile(data_path, 'eval');

   assert(exist(input_path, 'dir') == 7, ...
      'formal test input path does not exist: %s', input_path)
   assert(exist(eval_path, 'dir') == 7, ...
      'formal test eval path does not exist: %s', eval_path)

   if exist(output_path, 'dir') ~= 7
      mkdir(output_path);
   end

   % The formal suite expects these files under demo/data/input:
   %   met/met_kanm_kanm_2016_15m.mat
   %   met/met_kanl_kanl_2016_15m.mat
   %   spectral/kabs.mat
   %   spectral/kice.mat
   %   spectral/mie.mat
   %   spectral/solar.mat
   icemodel.config('ICEMODEL_DATA_PATH', data_path, ...
      'ICEMODEL_INPUT_PATH', input_path, ...
      'ICEMODEL_OUTPUT_PATH', output_path, ...
      'ICEMODEL_EVAL_PATH', eval_path, ...
      'ICEMODEL_USERDATA_PATH', fullfile(input_path, 'userdata'));
end
