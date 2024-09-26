function varargout = config(cfg)
   %CONFIG Configure icemodel project paths.
   %
   %  CFG = ICEMODEL.CONFIG()
   %  CFG = ICEMODEL.CONFIG(ICEMODEL_INPUT_PATH, pathname)
   %  CFG = ICEMODEL.CONFIG(ICEMODEL_OUTPUT_PATH, pathname)
   %
   % Description 
   %  
   %  CFG = ICEMODEL.CONFIG() Sets default icemodel environment variables
   %  including paths to model input files, specified by ICEMODEL_INPUT_PATH,
   %  and the location where model output files are saved, specified by
   %  ICEMODEL_OUTPUT_PATH.
   %  
   %  CFG = ICEMODEL.CONFIG(ICEMODEL_INPUT_PATH, pathname) Sets the environment
   %  variable ICEMODEL_INPUT_PATH, which specifies the path to the icemodel
   %  input files. Note that ICEMODEL_INPUT_PATH has a required sub-directory
   %  structure which is described below.
   %  
   %  CFG = ICEMODEL.CONFIG(ICEMODEL_OUTPUT_PATH, pathname) Sets the environment
   %  variable ICEMODEL_OUTPUT_PATH, which specifies the path where icemodel
   %  output files are saved.
   %
   %  CFG = ICEMODEL.CONFIG(ICEMODEL_CASE_NAME, casename) Sets a top-level
   %  icemodel project casename. If this value is provided, then the 
   %
   % Input arguments (name-value pairs)
   %
   %  ICEMODEL_CASE_NAME - an optional casename descriptor which can be used to
   %  organize simulation outputs. The default value is an empty string, which 
   %  ICEMODEL_INPUT_PATH - the path to the model input data. The default value
   %  is /path/to/this/repo/input.
   %  ICEMODEL_OUTPUT_PATH - the path where model output is saved. The default
   %  value is /path/to/this/repo/output.
   %  
   % Output arguments
   %  
   % ICEMODELDATAPATH - the path where model evaluation data is saved. The
   %  default location is ICEMODEL_INPUT_PATH/userdata.
   % 
   % The default model directory structure:
   %
   % icemodel/                            top level
   % icemodel/icemodel/                   model code
   % icemodel/input/                      input data (ICEMODEL_INPUT_PATH)
   % icemodel/input/met/                  input meteorological forcing data
   % icemodel/input/spectral/             input data for the spectral model
   % icemodel/input/userdata/             input user data (ICEMODELDATAPATH)
   % icemodel/output/                     model output data (ICEMODEL_OUTPUT_PATH)
   %
   % Note that .gitignore which ships with icemodel ignores input/ and output/.
   %
   % Examples
   %  
   %  1. Configure a default icemodel workspace.
   % 
   %  To configure a default icemodel workspace, call this function with no
   %  arguments:
   %     cfg = icemodel.config()
   %  
   %  This will set the environment variables ICEMODEL_INPUT_PATH and
   %  ICEMODEL_OUTPUT_PATH to the top-level `input` and `output` directories. If
   %  the `input` folder does not exist, the model will issue an error (input
   %  data is required to run the model). If the `output` directory does not
   %  exist, it will be created.
   %
   %  3. Specify custom input and output paths.
   %
   %     pathname_input = /path/to/input;
   %     pathname_output = /path/to/output
   %     config = icemodel.config("ICEMODEL_INPUT_PATH", pathname_input, ...
   %        "ICEMODEL_OUTPUT_PATH", pathname_output)
   % 
   % Note that in each case, 
   % 
   % See also: icemodel icemodel.setopts icemodel.run

   % Set default paths relative to the installation directory
   arguments
      cfg.ICEMODEL_CASE_NAME (1, :) = string.empty()
      cfg.ICEMODEL_INPUT_PATH (1, :) string = icemodel.internal.fullpath('input');
      cfg.ICEMODEL_OUTPUT_PATH (1, :) string = icemodel.internal.fullpath('output');
   end

   % Override the default options for case "demo". Note to users: This option is
   % included to create an isolated "demo/input" and "demo/output" directory
   % structure in the top-level icemodel path. 
   if cfg.ICEMODEL_CASE_NAME == "demo"
      demopath = icemodel.internal.fullpath("demo");
      cfg.ICEMODEL_INPUT_PATH = fullfile(demopath, "input");
      cfg.ICEMODEL_OUTPUT_PATH = fullfile(demopath, "output");
   end

   % Set userdata path relative to input/
   cfg.ICEMODELDATAPATH = fullfile(cfg.ICEMODEL_INPUT_PATH, 'userdata');

   % Toolbox metadata
   kwargs.ICEMODEL_VERSION = icemodel.internal.version();
   kwargs.ICEMODEL_REFERENCE = icemodel.internal.reference();
   kwargs.ICEMODEL_CONTACT = icemodel.internal.contact();

   % Set the environment variables
   for field = fieldnames(kwargs)'
      setenv(field{:}, kwargs.(field{:}))
   end

   % Return the config if requested
   if nargout == 1
      varargout{1} = kwargs;
   end
end
