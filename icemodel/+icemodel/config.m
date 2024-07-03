function varargout = config(varargin)
   %CONFIG Configure icemodel project paths.
   %
   %  CFG = ICEMODEL.CONFIG()
   %  CFG = ICEMODEL.CONFIG('ICEMODEL_INPUT_PATH', PATH_NAME)
   %  CFG = ICEMODEL.CONFIG('ICEMODEL_OUTPUT_PATH', PATH_NAME)
   %
   %% Description
   %
   %  CFG = ICEMODEL.CONFIG() Sets default icemodel environment variables,
   %  including paths to model input files (specified by ICEMODEL_INPUT_PATH),
   %  and the location where model output files are saved (specified by
   %  ICEMODEL_OUTPUT_PATH).
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
   %% Input arguments
   %
   %  ICEMODEL_INPUT_PATH  -  The path to the model input data. The default
   %                          value is /path/to/this/repo/input.
   %  ICEMODEL_OUTPUT_PATH -  The path where model output is saved. The default
   %                          value is /path/to/this/repo/output.
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
   %  ICEMODEL_INPUT_PATH  -  See description above.
   %  ICEMODEL_OUTPUT_PATH -  See description above.
   %  ICEMODEL_DATA_PATH   -  The path to additional (optional) user data files.
   %                          See the icemodel repo README for a description.
   %                          The path is set programmatically to
   %                          ICEMODEL_INPUT_PATH/userdata. Note that this
   %                          variable is primarily intended for user scripting
   %                          to conveniently access the data files saved in the
   %                          input/userdata folder. Users can place any files
   %                          inside userdata/ which may be of interest, and
   %                          programmatically access them using the
   %                          ICEMODEL_DATA_PATH environment variable.
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
   %  icemodel/input/                input data (ICEMODEL_INPUT_PATH)
   %  icemodel/input/met/            input meteorological forcing data
   %  icemodel/input/spectral/       input data for the spectral model
   %  icemodel/input/userdata/       input user data (ICEMODEL_DATA_PATH)
   %  icemodel/output/               model output data (ICEMODEL_OUTPUT_PATH)
   %
   % The default configuration above assumes ICEMODEL_INPUT_PATH and
   % ICEMODEL_OUTPUT_PATH are nested inside the top level project path,
   % alongside the model code, but this is not required. Set custom locations
   % for these folders by passing values for these variables to this function.
   %
   % Notes:
   %  - The .gitignore that ships with icemodel ignores input/ and output/.
   %  - ICEMODEL_INPUT_PATH can be set to any folder, but it must contain the
   %  following sub-directory structure, as outlined above:
   %
   %     input/met         (required)
   %     input/spectral    (required)
   %     input/userdata    (optional)
   %
   %  - See the example files in demo/input to understand how to create custom
   %  met data files, spectral data files, and alternative forcing data
   %  (userdata) files.
   %  - The ICEMODEL_OUTPUT_PATH can be set to any folder. When an individual
   %  simulation is configured using the icemodel.setopts() function, the
   %  following output subdirectory structure is created:
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
   %  This will set the environment variables ICEMODEL_INPUT_PATH and
   %  ICEMODEL_OUTPUT_PATH to the top-level `input` and `output` directories. If
   %  the `input` folder does not exist, the model will issue an error (input
   %  data is required to run the model). If the `output` directory does not
   %  exist, it will be created when a simulation is created with the
   %  icemodel.setopts function.
   %
   %  2. Specify custom input and output paths.
   %
   %     pathname_input = /path/to/input;
   %     pathname_output = /path/to/output;
   %     config = icemodel.config("ICEMODEL_INPUT_PATH", pathname_input, ...
   %        "ICEMODEL_OUTPUT_PATH", pathname_output)
   %
   % See also: icemodel icemodel.setopts icemodel.run

   % Input parsing
   kwargs = parseinputs(varargin{:});

   % Set userdata path relative to input/
   kwargs.ICEMODEL_DATA_PATH = fullfile(kwargs.ICEMODEL_INPUT_PATH, 'userdata');

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
   default_input_path = icemodel.internal.fullpath('input');
   default_output_path = icemodel.internal.fullpath('output');

   parser = inputParser();
   parser.addParameter('casename', char.empty(), @isscalartext)
   parser.addParameter('ICEMODEL_INPUT_PATH', default_input_path, @isscalartext)
   parser.addParameter('ICEMODEL_OUTPUT_PATH', default_output_path, @isscalartext)
   parser.parse(varargin{:})

   % Override the default options for case "demo". Note to users: This option is
   % included to create an isolated "demo/input" and "demo/output" directory
   % structure in the top-level icemodel path. The "casename" option is not used
   % for any other purpose.
   if parser.Results.casename == "demo"
      demopath = icemodel.internal.fullpath('demo');
      kwargs.ICEMODEL_INPUT_PATH = fullfile(demopath, 'input');
      kwargs.ICEMODEL_OUTPUT_PATH = fullfile(demopath, 'output');

   elseif ~isempty(parser.Results.casename)
      error('CASENAME argument currently only supports option "DEMO"')

   else
      kwargs.ICEMODEL_INPUT_PATH = parser.Results.ICEMODEL_INPUT_PATH;
      kwargs.ICEMODEL_OUTPUT_PATH = parser.Results.ICEMODEL_OUTPUT_PATH;
   end

   if ~(exist(kwargs.ICEMODEL_INPUT_PATH, 'dir') == 7)
      warning(['ICEMODEL_INPUT_PATH does not exist. ' ...
         'Create it and save the required input files there.'])
   end
end
