clean
%--------------------------------------------------------------------------
%   run dependent information
%--------------------------------------------------------------------------
savedata    =  true;
sitename    =  'region';
forcingdata =  'mar';
userdata    =  'modis';
uservars    =  'albedo';
meltmodel   =  'icemodel';

startyear   =  2008;
endyear     =  2018;
npts        =  1487;
si          =  671;              % 298, 671, 1046,  1420 
ei          =  744;             % 372, 744, 1116, 1487

% activate the right version
% activateicemodel('icemodel')

% MAKE SURE TO SET THE RIGHT PATH
drive          =  '/Users/mattcooper/mydata/';
opts.pathmet   =  [drive 'mar3.11/matfiles/region/level2/lists/'];
opts.pathsave  =  ['/Volumes/Samsung_T5b/icemodel/output/region/v10/' userdata '/'];
opts.pathuser  =  '/Users/mattcooper/myprojects/icemodel/input/';

% drive          =  '/Users/coop558/mydata/';
% opts.pathmet   =  [drive 'mar3.11/matfiles/region/level2/'];
% opts.pathsave  =  [drive 'icemodel/output/region/v10/skinmodel/' userdata '/'];
% opts.pathuser  =  '/Users/mattcooper/myprojects/icemodel/input/';

% to collapse region and version, need to use this format:
% opts.path.input    = setpath('GREENLAND/icemodel/input/');
% opts.path.output   = ['/Users/coop558/mydata/icemodel/output/v10/' ID '/'];
% opts.path.metdata  = [opts.path.input 'met/'];
% opts.path.initdata = [opts.path.input 'init/'];
% opts.path.userdata = [opts.path.input 'userData/'];


opts.savedata  =  savedata;

% drive          =  '/Users/coop558/mydata/';
% opts.pathmet   =  [drive 'mar3.11/matfiles/region/level2/' userdata '/'];
% opts.pathsave  =  [drive 'icemodel/model/experiment/v10/' userdata '/'];
% opts.pathuser  =  '/Users/coop558/MATLAB/GREENLAND/icemodel/input/';

% this is to use the new met files that have both mar and modis albedo
opts.userdata  =  userdata;

%--------------------------------------------------------------------------
%   run the model
%--------------------------------------------------------------------------

% initialize the model options
opts = a_opts_region(opts,sitename,meltmodel,startyear,endyear);

% display the run information
disp([meltmodel ' ' userdata ' ' int2str(startyear) ':'   ...
        int2str(endyear) ' ' int2str(si) ':' int2str(ei)])
     
% run the model for each point
for ipt = si:ei % 1045 = rio behar
    
   opts.ipt       = ipt;
   opts.metfname  = [opts.path.met 'met_' int2str(opts.ipt) '.mat'];
    
   % run the model
   [ice1,ice2] = icemodel_region(opts);

end

% save the opts
if ~exist([opts.pathsave 'opts/'],'dir'); mkdir([opts.pathsave 'opts/']); end
save([opts.pathsave 'opts/opts_' meltmodel '_' userdata '.mat'],'opts');
