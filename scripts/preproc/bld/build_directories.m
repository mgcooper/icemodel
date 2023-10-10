clean

% this is outdated but could be used as a template for setting up directories

% see /Users/coop558/mydata/icemodel/model/experiment/v8c/region/skinmodel
% need to create similar dir's for v10 icemodel for region runs

yri     =   2009;
yrf     =   2018;
nyears  =   yrf-yri+1;
basin   =   'ak4';

frc_met =   true; % force copy even if files already exist in destination
frc_fnc =   false; % functions 
frc_src =   false; % so
frc_drv =   false; % driver scripts
frc_spc =   false; % spectral data
frc_sol =   false; % solar data

% note: previously I had fmet = met1094.mat (or whatever the rio behar one
% is called) here, but I deleted and replaced with metYYYY in the loop

% the elements that need to be copied into each folder are 1) spectral
% profiles, 2) met files, 3) functions. The 

%------------------------------------------------------------------------------
% set paths to the folder where everything exists and the one you want to build
%------------------------------------------------------------------------------
p.src   = setpath('GREENLAND/icemodel/model/experiment/v8/behar/');
p.bld   = setpath(['GREENLAND/icemodel/model/experiment/v8/' basin '/']);
p.met   = setpath(['GREENLAND/icemodel/preproc/data/met/' basin '/']);
p.rad   = setpath('GREENLAND/icemodel/preproc/data/spectral/');

% NOTE: to separate the model from the location where the output is saved,             
%------------------------------------------------------------------------------
% make the directories
%------------------------------------------------------------------------------

if ~exist([p.bld],'dir')
    mkdir([p.bld]);
end

if ~exist([p.bld 'input'],'dir')
    mkdir([p.bld 'input']);
end

if ~exist([p.bld 'output'],'dir')
    mkdir([p.bld 'output']);
end

if ~exist([p.bld 'functions'],'dir')
    copyfile([p.src 'functions'],[p.bld 'functions']);
end

%------------------------------------------------------------------------------
% copy the static files
%------------------------------------------------------------------------------
flist   =   {   'a_icemodel_driver.m',          ...
                'b_icemodel_input.m',           ...
                'c_icemodel_initialize.m',      ...
                'd_icemodel.m',                 ...
                'e_icemodel_evaluation.m'       };

for n = 1:length(flist)
    f               =   flist{n};
    f_in            =   [p.src f];
    f_out           =   [p.bld f];
    
    if ~exist(f_out,'file') || frc_drv == true
        [status,msg,msgID] = copyfile(f_in,p.bld);
    end
end

%------------------------------------------------------------------------------
% copy the forcings
%------------------------------------------------------------------------------

for n = 1:nyears
    
    yy      =   int2str(yri+n-1);
    p.met_n =   [p.met yy '/'];
    p.out   =   [p.bld 'input/' yy '/'];
    fmet    =   ['met_' yy '.mat'];
    
  % build the directories (it's possible p.out exists, but not the others)
    if ~exist(p.out,'dir')
        mkdir(p.out);
    end
    if ~exist([p.out 'initialized/'],'dir')
        mkdir([p.out 'initialized/']);
    end
    if ~exist([p.out 'met/'],'dir')
        mkdir([p.out 'met/']);
    end
    if ~exist([p.out 'spectral/'],'dir')
        mkdir([p.out 'spectral/']);
    end
    
  % copy the met files
    if ~exist([p.out 'met/' fmet],'file') || frc_met == true
        copyfile([p.met_n fmet],[p.out 'met/']);
    end
    
  % copy the spectral profiles
    if ~exist([p.out 'spectral/kabs.mat'],'file') || frc_spc == true
        copyfile([p.rad 'kabs.mat'],[p.out 'spectral/']);
    end
    
    if ~exist([p.out 'spectral/kice.mat'],'file') || frc_spc == true
        copyfile([p.rad 'kice.mat'],[p.out 'spectral/']);
    end
    
    if ~exist([p.out 'spectral/mie.mat'],'file') || frc_spc == true
        copyfile([p.rad 'mie.mat'],[p.out 'spectral/']);
    end
    
    if ~exist([p.out 'spectral/solar.mat'],'file') || frc_sol == true
        copyfile([p.rad 'solar.mat'],[p.out 'spectral/']);
    end
    
end

%------------------------------------------------------------------------------
% make an output folder for each year
%------------------------------------------------------------------------------
for n = 1:nyears
    
    yy = int2str(yri+n-1);
    
    if ~exist([p.bld 'output/' yy],'dir')
        mkdir([p.bld 'output/' yy]);
    end
end
    
    