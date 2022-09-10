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

startyear   =   2008;
endyear     =   2018;
npts        =   1487;
si          =   1;              % 298, 671, 1046,  1420 
ei          =   743;             % 372, 744, 1116, 1487

% warning off
% rmpath(genpath('/Volumes/GoogleDrive/My Drive/MATLAB/GREENLAND/icemodel/model'));
% modelpath = '/Users/mattcooper/myprojects/icemodel/skinmodel';
% addpath(genpath(modelpath));
% cd(modelpath)
% warning on

% MAKE SURE TO SET THE RIGHT PATH
drive             =  '/Users/mattcooper/mydata/';
opts.path.met     =  [drive 'mar3.11/matfiles/region/level2/lists/'];
opts.path.input   =  '/Users/mattcooper/myprojects/icemodel/input/';
opts.path.output  =  '/Volumes/Samsung_T5b/icemodel/output/region/v10/';
opts.path.output  =  [opts.path.output meltmodel '/' userdata '/'];

% drive             =  '/Users/coop558/mydata/';
% opts.path.met     =  [drive 'mar3.11/matfiles/region/level2/'];
% opts.path.input   =  setpath('GREENLAND/icemodel/input/');

opts.savedata  =  savedata;

% this is to use the new met files that have both mar and modis albedo
opts.userdata  =  userdata;

%% run the model

disp([meltmodel ' ' userdata ' ' int2str(startyear) ':'   ...
        int2str(endyear) ' ' int2str(si) ':' int2str(ei)])

% initialize the model options
opts = a_opts_region(opts,sitename,meltmodel,startyear,endyear);

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
