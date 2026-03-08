function varargout = config(varargin)
   %CONFIG Configure icemodel project paths.
   %
   %  CFG = ICEMODEL.CONFIG()
   %  CFG = ICEMODEL.CONFIG('ICEMODEL_DATA_PATH', PATH_NAME)
   %  CFG = ICEMODEL.CONFIG('ICEMODEL_INPUT_PATH', PATH_NAME)
   %  CFG = ICEMODEL.CONFIG('ICEMODEL_OUTPUT_PATH', PATH_NAME)
   %  CFG = ICEMODEL.CONFIG('ICEMODEL_EVAL_PATH', PATH_NAME)
   %  CFG = ICEMODEL.CONFIG('ICEMODEL_USERDATA_PATH', PATH_NAME)
   %
   %% Description
   %
   %  CFG = ICEMODEL.CONFIG() Sets default icemodel environment variables,
   %  including the top-level data root (ICEMODEL_DATA_PATH), paths to model
   %  input files (ICEMODEL_INPUT_PATH), evaluation/reference data
   %  (ICEMODEL_EVAL_PATH), optional userdata swap files
   %  (ICEMODEL_USERDATA_PATH), and the location where model output files are
   %  saved (ICEMODEL_OUTPUT_PATH).
   %
   %  CFG = ICEMODEL.CONFIG("ICEMODEL_DATA_PATH", PATH_NAME) Sets the
   %  environment variable ICEMODEL_DATA_PATH to PATH_NAME. PATH_NAME is a
   %  scalar text object (string or row vector of characters) which specifies
   %  the full path where the default icemodel data directory structure exists.
   %  Unless overridden individually, the other path variables are derived from
   %  this root as:
   %
   %     ICEMODEL_INPUT_PATH    = ICEMODEL_DATA_PATH/input
   %     ICEMODEL_OUTPUT_PATH   = ICEMODEL_DATA_PATH/output
   %     ICEMODEL_EVAL_PATH     = ICEMODEL_DATA_PATH/eval
   %     ICEMODEL_USERDATA_PATH = ICEMODEL_INPUT_PATH/userdata
   %
   %  CFG = ICEMODEL.CONFIG("ICEMODEL_INPUT_PATH", PATH_NAME) Sets the
   %  environment variable ICEMODEL_INPUT_PATH to PATH_NAME. PATH_NAME is a
   %  scalar text object (string or row vector of characters) which specifies
   %  the full path where the icemodel input files exist. Note that the
   %  ICEMODEL_INPUT_PATH sub-directory structure must conform to a specific
   %  format described below in the "Model directory structure" section.
   %
   %  CFG = ICEMODEL.CONFIG("ICEMODEL_OUTPUT_PATH", PATH_NAME) Sets the
   %  environment variable ICEMODEL_OUTPUT_PATH to PATH_NAME. PATH_NAME is a
   %  scalar text object (string or row vector of characters) which specifies
   %  the full path where the icemodel output files are saved. Note that the
   %  ICEMODEL_OUTPUT_PATH sub-directory structure is created automatically when
   %  running an icemodel simulation, described below in the "Model directory
   %  structure" section.
   %
   %  CFG = ICEMODEL.CONFIG("ICEMODEL_EVAL_PATH", PATH_NAME) Sets the
   %  environment variable ICEMODEL_EVAL_PATH to PATH_NAME. PATH_NAME is a
   %  scalar text object (string or row vector of characters) which specifies
   %  the full path where evaluation/reference data are stored, such as
   %  ablation observations or other external comparison datasets.
   %
   %  CFG = ICEMODEL.CONFIG("ICEMODEL_USERDATA_PATH", PATH_NAME) Sets the
   %  environment variable ICEMODEL_USERDATA_PATH to PATH_NAME. PATH_NAME is a
   %  scalar text object (string or row vector of characters) which specifies
   %  the full path where optional userdata swap files are stored.
   %
   %% Input arguments
   %
   %  ICEMODEL_DATA_PATH     - The top-level data path. The default value is
   %                           /path/to/this/repo/data.
   %  ICEMODEL_INPUT_PATH    - The path to the model input data. The default
   %                           value is ICEMODEL_DATA_PATH/input.
   %  ICEMODEL_OUTPUT_PATH   - The path where model output is saved. The
   %                           default value is ICEMODEL_DATA_PATH/output.
   %  ICEMODEL_EVAL_PATH     - The path to evaluation/reference data. The
   %                           default value is ICEMODEL_DATA_PATH/eval.
   %  ICEMODEL_USERDATA_PATH - The path to optional userdata swap files. The
   %                           default value is ICEMODEL_INPUT_PATH/userdata.
   %
   % Note that all input arguments are specified as name-value pairs, also known
   % as keyword arguments (kwargs). Specify inputs using comma-separated
   % name-value pairs: icemodel.config("ICEMODEL_INPUT_PATH", PATH_NAME), or
   % Name=Value syntax: icemodel.config(ICEMODEL_INPUT_PATH=PATH_NAME).
   %
   %% Output arguments
   %
   %  CFG - A structure containing the configuration variable name-value pairs.
   %
   %  CFG contains the following name-value pairs which specify icemodel paths:
   %  ICEMODEL_DATA_PATH     - See description above.
   %  ICEMODEL_INPUT_PATH    - See description above.
   %  ICEMODEL_OUTPUT_PATH   - See description above.
   %  ICEMODEL_EVAL_PATH     - See description above.
   %  ICEMODEL_USERDATA_PATH - See description above.
   %
   %  CFG also contains the following toolbox metadata variables which are set
   %  automatically by internal functions. These are not settable by users but
   %  can be accessed for reference.
   %
   %  ICEMODEL_VERSION     -  The model version number. Set automatically by
   %                          icemodel.internal.version().
   %  ICEMODEL_REFERENCE   -  The icemodel code reference for citations. Set
   %                          automatically by icemodel.internal.reference().
   %  ICEMODEL_CONTACT     -  Contact information for questions about icemodel.
   %                          Set automatically by icemodel.internal.contact()
   %
   %% Model directory structure
   %
   % The following directory structure demonstrates the default configuration:
   %
   %  icemodel/                      top level project path (this repo)
   %  icemodel/icemodel/             model code
   %  icemodel/data/                 top-level data root (ICEMODEL_DATA_PATH)
   %  icemodel/data/input/           input data (ICEMODEL_INPUT_PATH)
   %  icemodel/data/input/met/       input meteorological forcing data
   %  icemodel/data/input/spectral/  input data for the spectral model
   %  icemodel/data/input/userdata/  input user data (default ICEMODEL_USERDATA_PATH)
   %  icemodel/data/eval/            evaluation/reference data (ICEMODEL_EVAL_PATH)
   %  icemodel/data/output/          model output data (ICEMODEL_OUTPUT_PATH)
   %
   % The default configuration above assumes ICEMODEL_DATA_PATH is nested inside
   % the top level project path, alongside the model code, but this is not
   % required. Set custom locations for any of these folders by passing values
   % for the corresponding variables to this function.
   %
   % Notes:
   %  - The .gitignore that ships with icemodel ignores data/input/,
   %    data/eval/, and data/output/.
   %  - ICEMODEL_INPUT_PATH can be set to any folder, but it must contain the
   %  following sub-directory structure, as outlined above:
   %
   %     input/met         (required)
   %     input/spectral    (required)
   %
   %  - ICEMODEL_USERDATA_PATH can be set independently, but by default it is
   %    assumed to be ICEMODEL_INPUT_PATH/userdata.
   %
   %  - See the example files in demo/data/input to understand how to create
   %    custom met data files, spectral data files, and alternative forcing data
   %    (userdata) files.
   %  - The ICEMODEL_OUTPUT_PATH can be set to any folder. When an individual
   %    simulation is run, the following output directory structure is created:
   %
   %     ICEMODEL_OUTPUT_PATH/<sitename>/<smbmodel>
   %
   %  Here, <sitename> and <smbmodel> are variables passed to icemodel.setopts,
   %  which specify the site name for the simulation and which surface mass
   %  balance model (smbmodel) to run: SkinModel or IceModel.
   %
   %  It is often desirable to add an additional sub-folder (e.g., to set up
   %  different runs for the same sitename and smbmodel). To do so, pass
   %  the optional <testname> argument to icemodel.setopts, and the following
   %  output subdirectory structure will be created:
   %
   %     ICEMODEL_OUTPUT_PATH/<sitename>/<smbmodel>/<testname>
   %
   %  See demo.m for an example of how to set simulation options using
   %  icemodel.setopts.
   %
   %% Examples
   %
   %  1. Configure a default icemodel workspace.
   %
   %  To configure a default icemodel workspace, call this function with no
   %  arguments:
   %     cfg = icemodel.config()
   %
   %  This will set the environment variables ICEMODEL_DATA_PATH,
   %  ICEMODEL_INPUT_PATH, ICEMODEL_EVAL_PATH, ICEMODEL_USERDATA_PATH, and
   %  ICEMODEL_OUTPUT_PATH to the default directories under `data/`. If the
   %  `input` folder does not exist, the model will issue an error (input data
   %  is required to run the model). If the `output` directory does not exist,
   %  it will be created automatically when a simulation is run.
   %
   %  2. Specify custom input and output paths.
   %
   %     pathname_data = /path/to/data;
   %     pathname_output = /path/to/output;
   %     config = icemodel.config("ICEMODEL_DATA_PATH", pathname_data, ...
   %        "ICEMODEL_OUTPUT_PATH", pathname_output)
   %
   % See also: icemodel icemodel.setopts icemodel.run.point

   % Input parsing
   kwargs = parseinputs(varargin{:});

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

function kwargs = parseinputs(varargin)

   % Set default paths relative to the installation directory
   default_data_path = icemodel.internal.fullpath('data');
   default_input_path = fullfile(default_data_path, 'input');
   default_output_path = fullfile(default_data_path, 'output');
   default_eval_path = fullfile(default_data_path, 'eval');
   default_user_path = fullfile(default_input_path, 'userdata');

   parser = inputParser();
   parser.addParameter('casename', char.empty(), @isscalartext)
   parser.addParameter('ICEMODEL_DATA_PATH', default_data_path, @isscalartext)
   parser.addParameter('ICEMODEL_INPUT_PATH', default_input_path, @isscalartext)
   parser.addParameter('ICEMODEL_OUTPUT_PATH', default_output_path, @isscalartext)
   parser.addParameter('ICEMODEL_EVAL_PATH', default_eval_path, @isscalartext)
   parser.addParameter('ICEMODEL_USERDATA_PATH', default_user_path, @isscalartext)
   parser.parse(varargin{:})

   % Override the default options for case "demo". Note to users: This option
   % is used to create an isolated "demo/data/..." directory structure in the
   % top-level icemodel path. There is only one "casename" option - "demo" -
   % its only purpose is to configure the demo.
   if strcmp(parser.Results.casename, 'demo')
      default_data_path = icemodel.internal.fullpath('demo', 'data');
      default_input_path = fullfile(default_data_path, 'input');
      default_output_path = fullfile(default_data_path, 'output');
      default_eval_path = fullfile(default_data_path, 'eval');
      default_user_path = fullfile(default_input_path, 'userdata');
   elseif ~isempty(parser.Results.casename)
      error('CASENAME argument currently only supports option "DEMO"')
   end

   kwargs.ICEMODEL_DATA_PATH = getresult(parser, 'ICEMODEL_DATA_PATH', ...
      default_data_path);
   kwargs.ICEMODEL_INPUT_PATH = getresult(parser, 'ICEMODEL_INPUT_PATH', ...
      default_input_path);
   kwargs.ICEMODEL_OUTPUT_PATH = getresult(parser, 'ICEMODEL_OUTPUT_PATH', ...
      default_output_path);
   kwargs.ICEMODEL_EVAL_PATH = getresult(parser, 'ICEMODEL_EVAL_PATH', ...
      default_eval_path);
   kwargs.ICEMODEL_USERDATA_PATH = getresult(parser, 'ICEMODEL_USERDATA_PATH', ...
      default_user_path);

   if ~(exist(kwargs.ICEMODEL_INPUT_PATH, 'dir') == 7)
      warning(['ICEMODEL_INPUT_PATH does not exist. ' ...
         'Create it and save the required input files there.'])
   end
end

function value = getresult(parser, name, default_value)
   if ismember(name, parser.UsingDefaults)
      value = default_value;
   else
      value = parser.Results.(name);
   end
end
