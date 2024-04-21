function out = config()
   %CONFIG Configure project preferences
   %
   %  CFG = ICEMODEL.CONFIG() Configures project preferences including paths to
   %  data files and the location of model output files.
   %
   % Description
   %
   %  This is an example configuration. The model needs to know where the input
   %  data is located, and (optionally) where to save the output data. this
   %  script shows one way to accomplish that. more information is provided
   %  below.
   %
   % Proerties
   %
   %  ICEMODELPATH - the top-level model path (see 'more information' below)
   %  ICEMODELDATAPATH - the path to the evaluation data
   %  ICEMODELINPUTPATH - the path to the input met data
   %  ICEMODELOUTPUTPATH - the path where model output is saved
   %  ICEMODELUSERPATH - an additional path that is used for accessing large
   %  datasets stored on external hard drives, or anything else the user wishes
   %  to use it for. It is not used in the model, but may be used in some
   %  plotting scripts or pre-processing scripts.
   %
   %  NOTE: ICEMODELPATH is not currently used, and ICEMODELOUTPUTPATH is only
   %  used in the icemodel_run script. ICEMODELINPUTPATH should point to input/
   %  and ICEMODELDATAPATH should point to input/userdata )
   %
   % See also: icemodel icemodel.setopts icemodel.run

   HOMEPATH = getenv('HOME');

   cfg.ICEMODELPATH = fullfile(HOMEPATH, 'myprojects/matlab/icemodel');
   cfg.ICEMODELINPUTPATH = fullfile(HOMEPATH,'myprojects/matlab/runoff/data/icemodel/input');
   cfg.ICEMODELDATAPATH = fullfile(HOMEPATH, 'myprojects/matlab/runoff/data/icemodel/eval');
   cfg.ICEMODELOUTPUTPATH = '/Volumes/Samsung_T5b/icemodel/output/v10b';
   cfg.ICEMODELUSERPATH = '/Volumes/Samsung_T5b';

   cfg.ICEMODEL_VERSION = icemodel.internal.version();
   cfg.ICEMODEL_REFERENCE = icemodel.internal.reference();
   cfg.ICEMODEL_CONTACT = icemodel.internal.contact();

   for field = fieldnames(cfg)'
      setenv(field{:}, cfg.(field{:}))
   end

   if nargout > 0
      out = cfg;
   end
end

% ----------------
% MORE INFORMATION
% ----------------

% this is the default model directory structure:

% icemodel/                            top level
% icemodel/drive/                      scripts to configure and run the model
% icemodel/functions/                  model code
% icemodel/input/                      input data
% icemodel/input/met/                  input meteorological forcing data
% icemodel/input/spectral/             input spectrally-dependent data
% icemodel/input/userdata/             input user data
% icemodel/output/                     output data (produced by the model)
% icemodel/scripts/                    additional scripts

% (note: the .gitignore that ships with icemodel ignores input/ and output/)

% as mentioned, the model needs to know where the input forcing data is stored.
% one option is to run the model from the top-level folder, in which case the
% model will look for the met data forcing file in the input/met folder.

% a second option is to define environment variables that define the path to
% the input data and/or output data, as above. if these environment variables
% exist, the model will look for the input data in ICEMODELINPUTPATH first and
% if not found, it will look for the data in $PWD/input/met/. Similarly, if
% ICEMODELOUTPUTPATH is defined, it will save the data there, otherwise it will
% look for $PWD/output/. if that fails, it will create a temporary directory
% $PWD/icemodel_output_tmp/ and save the data there.

% a third option is to simply use the 'addpath' function to ensure the input
% data exists somewhere on your path, and the model will look for it and, if
% found, will build the input path accordingly. In this case, the 'input/'
% folder must also contain the 'met/', 'spectral/', and 'userdata/'
% subdirectories as indicated in the directory structure above.

% the reason you might choose one option or another depends on how you use the
% model. in most cases, it should make sense to use the default directory
% structure and/or define the environment variables above, either here, or at
% the top of your drive script. you might also define them in startup.m, and
% override them here or in your drive script as needed, for example if you run
% lots of simulations that involve hundreds or more met files and generate many
% GBs of data, you may want to save the input and output data somewhere other
% than where the model code is saved (as i have it setup above).


% setenv('ICEMODELDATAPATH',fullfile(HOMEPATH,'myprojects/matlab/runoff/data/icemodel/eval'));
% setenv('ICEMODELINPUTPATH',fullfile(HOMEPATH,'myprojects/matlab/runoff/data/icemodel/input'));