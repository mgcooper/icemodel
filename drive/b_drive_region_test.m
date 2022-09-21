clean
%--------------------------------------------------------------------------
%   run dependent information
%--------------------------------------------------------------------------
workoff skinmodel
% workon runoff
workon icemodel


savedata    =  true;
sitename    =  'region';
forcingdata =  'mar';
userdata    =  'mar';
uservars    =  'albedo';
meltmodel   =  'icemodel';

startyear   =   2008;
endyear     =   2018;

% load the bad rh test data to find the worst offenders
load(setpath('runoff/NEW/badrh_2.mat','project')); badrh = test; clear test;
% ibad = find(badrh.numbad==max(badrh.numbad));
ibadold  = find(badrh.numover>4);         % indices of the old metfiles 
ibadnew  = findnewidxfromold(ibadold);    % indices of the new metfiles
npts     = numel(ibadold);

% MAKE SURE TO SET THE RIGHT PATH
drive             =  '/Volumes/Samsung_T5b/';
opts.path.met     =  setpath('mar3.11/matfiles/region/level2/lists3/',drive);
opts.path.input   =  setpath('runoff/data/icemodel/input/','project');
opts.path.output  =  setpath('icemodel/output/region/v10_b/',drive);
opts.path.output  =  [opts.path.output meltmodel '/' userdata '/'];

opts.savedata     =  savedata;
opts.userdata     =  userdata; % to set mar/modis albedo


%% run the model

% disp([meltmodel ' ' userdata ' ' int2str(startyear) ':'   ...
%         int2str(endyear) ' ' int2str(si) ':' int2str(ei)])

% initialize the model options
opts = a_opts_region(opts,sitename,meltmodel,startyear,endyear);

% run the model for each point
for n = 1:numel(ibadnew)

   ipt            = ibadnew(n);
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