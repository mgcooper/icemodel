function [ice1, ice2, met, opts] = point(kwargs)
   %POINT Run one point-scale simulation and post-process its output.
   %
   %  [ice1, ice2, met, opts] = icemodel.run.point()
   %  [ice1, ice2, met, opts] = icemodel.run.point("sitename", "kanm", ...
   %     "simyears", 2016, "smbmodel", "icemodel")
   %  [ice1, ice2, met, opts] = icemodel.run.point(_, ...
   %     "overrides", struct('dt', 900, 'solver', 3))
   %
   %  The arguments below select the case. OVERRIDES holds any other option
   %  by field name, and icemodel.resetopts applies it to the options this
   %  function builds, so a caller sets a non-standard value without editing
   %  icemodel.setopts. An override wins over the case arguments.
   %
   %  ice1 and ice2 are the post-processed model output, met is the forcing
   %  timetable, and opts is the resolved options struct of the run.
   %
   % See also: icemodel.setopts icemodel.resetopts icemodel.configureRun

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
      kwargs.overrides (1, 1) struct = struct()
   end
   [saveflag, sitename, forcings, userdata, uservars, ...
      smbmodel, simyears, gridcell, testname, backupflag, n_spinup_years, ...
      overrides] = deal(kwargs.saveflag, kwargs.sitename, kwargs.forcings, ...
      kwargs.userdata, kwargs.uservars, kwargs.smbmodel, kwargs.simyears, ...
      kwargs.gridcell, kwargs.testname, kwargs.backupflag, ...
      kwargs.n_spinup_years, kwargs.overrides);

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

   % Apply the caller overrides last, so an override wins over a case
   % argument and over the gridcell metfname. resetopts clears the derived
   % fields that depend on each override, and configureRun rebuilds them
   % inside the model entry point.
   if ~isempty(fieldnames(overrides))
      override_args = namedargs2cell(overrides);
      opts = icemodel.resetopts(opts, override_args{:});
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
      [ice1, ice2, met] = icemodel.postprocess( ...
         ice1, ice2, opts, opts.output_years);
   end
end
