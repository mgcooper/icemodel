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

% MAKE SURE TO SET THE RIGHT PATH
drive          =  '/Users/mattcooper/data/';
% opts.pathmet   =  [drive 'mar3.11/matfiles/region/level2/old/'];
opts.pathmet   =  [drive 'mar3.11/matfiles/region/level2/lists/'];
opts.pathsave  =  ['/Volumes/Samsung_T5b/icemodel/output/region/v10/' userdata '/'];
opts.pathuser  =  '/Users/mattcooper/projects/icemodel/input/';

% drive          =  '/Users/coop558/mydata/';
% opts.pathmet   =  [drive 'mar3.11/matfiles/region/level2/' userdata '/'];
% opts.pathsave  =  [drive 'icemodel/model/experiment/v10/' userdata '/'];
% opts.pathuser  =  '/Users/coop558/MATLAB/GREENLAND/icemodel/input/';

opts.savedata  =  savedata;

% this is unique to this version, which uses the new met files that have 
% both mar and modis albedo
opts.userdata  =  userdata;

%% run the model

disp([meltmodel ' ' userdata ' ' int2str(startyear) ':'   ...
        int2str(endyear) ' ' int2str(si) ':' int2str(ei)])

% initialize the model options
opts = a_opts_region(opts,sitename,meltmodel,startyear,endyear);

% run the model for each point
for ipt = si:ei % 1045 = rio behar
    
   opts.ipt    = ipt;
    
   % run the model
   [ice1,ice2] = icemodel_region(opts);

end

% save the opts
if ~exist([opts.pathsave 'opts/'],'dir'); mkdir([opts.pathsave 'opts/']); end
save([opts.pathsave 'opts/opts_' meltmodel '_' userdata '.mat'],'opts');

