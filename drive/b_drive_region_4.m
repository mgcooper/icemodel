clean
%--------------------------------------------------------------------------
%   run dependent information
%--------------------------------------------------------------------------
workoff skinmodel
workon icemodel

savedata    =  true;
sitename    =  'region';
forcingdata =  'mar';
userdata    =  'modis';
uservars    =  'albedo';
meltmodel   =  'icemodel';

startyear   =   2008;
endyear     =   2018;

% npts        =   1487;
% si          =   1;              % 298, 671, 1046,  1420 
% ei          =   743;             % 372, 744, 1116, 1487

load('private/idxnew.mat','idxnew')
npts        =   numel(idxnew);
si          =   497;              % 298, 671, 1046,  1420 
ei          =   npts-196;             % 372, 744, 1116, 1487
% ei          =   npts;             % 372, 744, 1116, 1487

% npts        =  2479;
% si          =  1;
% ei          =  1239;

% MAKE SURE TO SET THE RIGHT PATH
drive             =  '/Volumes/Samsung_T5b/';
opts.path.met     =  setpath('mar3.11/matfiles/region/level2/lists3/',drive);
opts.path.input   =  setpath('runoff/data/icemodel/input/','project');
opts.path.output  =  setpath('icemodel/output/region/v10_b/',drive);
opts.path.output  =  [opts.path.output meltmodel '/' userdata '/'];

opts.savedata     =  savedata;
opts.userdata     =  userdata; % to set mar/modis albedo

%% run the model

disp([meltmodel ' ' userdata ' ' int2str(startyear) ':'   ...
        int2str(endyear) ' ' int2str(si) ':' int2str(ei)])

% initialize the model options
opts = a_opts_region(opts,sitename,meltmodel,startyear,endyear);

% run the model for each point
% for ipt = si:ei % 1045 = rio behar
for n = si:ei

   ipt            = idxnew(n);
   opts.ipt       = ipt;
   opts.metfname  = [opts.path.met 'met_' int2str(opts.ipt) '.mat'];
    
   % run the model
   [ice1,ice2] = icemodel_region(opts);

end

% save the opts
if ~exist([opts.path.output 'opts/'],'dir')
   mkdir([opts.path.output 'opts/']); 
end
save([opts.path.output 'opts/opts_' meltmodel '_' userdata '.mat'],'opts');