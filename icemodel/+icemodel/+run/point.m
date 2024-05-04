function [ice1, ice2, met, opts] = point(kwargs)

   arguments (Input)
      kwargs.saveflag (1, 1) logical = false
      kwargs.sitename (1, :) string {mustBeMember(kwargs.sitename, ...
         ["sector", "behar", "ak4", "slv1", "slv2", "upperbasin", ...
         "kanm", "kanl"])} = "kanm"
      kwargs.forcings (1, :) string {mustBeMember(kwargs.forcings, ...
         ["mar", "kanm", "kanl"])} = "kanm"
      kwargs.userdata (1, :) string {mustBeMember(kwargs.userdata, ...
         ["mar", "modis", "merra", "racmo", "kanm", "kanl"])} = []
      kwargs.uservars (1, :) string {mustBeMember(kwargs.uservars, ...
         ["albedo", "tair", "swd", "lwd", "rh", "wsdp"])} = []
      kwargs.smbmodel (1, :) string {mustBeMember(kwargs.smbmodel, ...
         ["icemodel", "skinmodel"])} = "skinmodel"
      kwargs.simyears (1, :) double = []
      kwargs.gridcell (1, :) double = []
      kwargs.testname (1, :) string = []
      kwargs.backupflag (1, 1) logical = false
   end
   [saveflag, sitename, forcings, userdata, uservars, ...
      smbmodel, simyears, gridcell, testname, backupflag] ...
      = deal(kwargs.saveflag, kwargs.sitename, kwargs.forcings, ...
      kwargs.userdata, kwargs.uservars, kwargs.smbmodel, kwargs.simyears, ...
      kwargs.gridcell, kwargs.testname, kwargs.backupflag);

   if isempty(userdata)
      userdata = forcings;
   end

   % varargin = struct2cell(kwargs);
   % [varargin{:}] = convertStringsToChars(varargin{:});

   %% Set the model options
   opts = icemodel.setopts(smbmodel, sitename, simyears, forcings, ...
      userdata, uservars, saveflag, testname, backupflag);

   if notempty(gridcell)
      opts.metfname = {fullfile(opts.pathinput, 'met', 'sector', ...
         ['met_' int2str(gridcell) '.mat'])};
   end

   % run the model
   switch smbmodel
      case 'icemodel'
         tic; [ice1, ice2] = icemodel(opts); toc
      case 'skinmodel'
         tic; [ice1, ice2] = skinmodel(opts); toc
   end

   % load the met data and run the post processing
   if saveflag
      [ice1, ice2, met] = icemodel.loadresults(opts);
   else
      [ice1, ice2, met] = POSTPROC(ice1, ice2, opts, simyears);
   end
end
