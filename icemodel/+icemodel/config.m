function varargout = config(kwargs, cfg)
   %CONFIG Configure project preferences
   %
   %  CFG = ICEMODEL.CONFIG() Configures project preferences including paths to
   %  data files and the location of model output files.
   %
   % Description
   %
   %  This is an example configuration. The model needs to know where the input
   %  data is located, and (optionally) where to save the output data. This
   %  script shows one way to accomplish that.
   %
   % Properties
   %
   %  ICEMODELPATH - the top-level model path
   %  ICEMODELINPUTPATH - the path to the model input data
   %  ICEMODELOUTPUTPATH - the path where model output is saved
   %  ICEMODELUSERPATH - an additional path used for large datasets (e.g.,
   %  external hard drive), or anything else.
   %
   % See also: icemodel icemodel.setopts icemodel.run

   arguments
      kwargs.casename (1, :) = []
      cfg.ICEMODELPATH (1, :) string = icemodel.internal.fullpath();
      cfg.ICEMODELINPUTPATH (1, :) string = icemodel.internal.fullpath('input');
      cfg.ICEMODELOUTPUTPATH (1, :) string = icemodel.internal.fullpath('output');
   end

   % Override the default options for case "demo"
   if kwargs.casename == "demo"
      demopath = icemodel.internal.fullpath("demo");
      cfg.ICEMODELINPUTPATH = fullfile(demopath, "input");
      cfg.ICEMODELOUTPUTPATH = fullfile(demopath, "output");
   end

   cfg.ICEMODELDATAPATH = fullfile(cfg.ICEMODELINPUTPATH, 'userdata');
   cfg.ICEMODEL_VERSION = icemodel.internal.version();
   cfg.ICEMODEL_REFERENCE = icemodel.internal.reference();
   cfg.ICEMODEL_CONTACT = icemodel.internal.contact();

   for field = fieldnames(cfg)'
      setenv(field{:}, cfg.(field{:}))
   end

   if nargout == 1
      varargout{1} = cfg;
   end
end

% The default model directory structure:

% icemodel/                            top level
% icemodel/icemodel/                   model code
% icemodel/input/                      input data
% icemodel/input/met/                  input meteorological forcing data
% icemodel/input/spectral/             input data for the spectral model
% icemodel/input/userdata/             input user data
% icemodel/output/                     output data (produced by the model)

% note that the .gitignore that ships with icemodel ignores input/ and output/
