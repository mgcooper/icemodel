function [ice1, ice2, met, opts] = point(kwargs)

   arguments (Input)
      kwargs.saveflag (1, 1) logical = false
      kwargs.sitename (1, :) string { ...
         icemodel.validators.mustBeSiteName(kwargs.sitename)} = "kanm"
      kwargs.forcings (1, :) string { ...
         icemodel.validators.mustBeForcingName(kwargs.forcings)} = "kanm"
      kwargs.userdata (1, :) string { ...
         icemodel.validators.mustBeUserdataName(kwargs.userdata)} = []
      kwargs.uservars (1, :) string { ...
         icemodel.validators.mustBeUservarName(kwargs.uservars)} = []
      kwargs.smbmodel (1, :) string { ...
         icemodel.validators.mustBeSmbmodelName(kwargs.smbmodel)} = "skinmodel"
      kwargs.simyears (1, :) double = []
      kwargs.gridcell (1, :) double = []
      kwargs.testname (1, :) string = []
      kwargs.backupflag (1, 1) logical = false
      kwargs.n_spinup_years (1, 1) double {mustBeNonnegative, mustBeInteger} = 0
   end
   [saveflag, sitename, forcings, userdata, uservars, ...
      smbmodel, simyears, gridcell, testname, backupflag, n_spinup_years] ...
      = deal(kwargs.saveflag, kwargs.sitename, kwargs.forcings, ...
      kwargs.userdata, kwargs.uservars, kwargs.smbmodel, kwargs.simyears, ...
      kwargs.gridcell, kwargs.testname, kwargs.backupflag, kwargs.n_spinup_years);

   if isempty(userdata)
      userdata = forcings;
   end

   % varargin = struct2cell(kwargs);
   % [varargin{:}] = convertStringsToChars(varargin{:});

   %% Set the model options
   opts = icemodel.setopts(smbmodel, sitename, simyears, forcings, ...
      userdata, uservars, testname, saveflag, backupflag, ...
      'n_spinup_years', n_spinup_years);

   if notempty(gridcell)
      opts = icemodel.resetopts(opts, 'metfname', ...
         {fullfile(opts.pathinput, 'met', 'sector', ...
         ['met_' int2str(gridcell) '.mat'])});
   end

   % run the model
   switch smbmodel
      case 'icemodel'
         tic; [ice1, ice2, opts] = icemodel(opts); toc
      case 'skinmodel'
         tic; [ice1, ice2, opts] = skinmodel(opts); toc
   end

   % load the met data and run the post processing
   if saveflag
      [ice1, ice2, met] = icemodel.loadresults(opts);
   else
      [ice1, ice2, met] = POSTPROC(ice1, ice2, opts, opts.output_years);
   end
end
