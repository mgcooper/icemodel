function opts = icemodel_opts(opts,sitename,simyear,meltmodel,forcingdata,  ...
   userdata,uservars,startyear,endyear)

opts.sitename        =  sitename;
opts.yyyy            =  simyear;

% model options
opts.annual_loops    =  1;
opts.use_init        =  false;
opts.calendar_type   =  'noleap';
opts.dt              =  900.0;  % timestep                          [s]
opts.kthermal        =  12;     % thermal conductivity
opts.z_0             =  0.001;  % Surface aero. roughness length    [m]
opts.z_obs           =  3.0;    % Height of wind and air temp obs   [m]
opts.kabs_user       =  true;
opts.tlagcolumn      =  6*3600/opts.dt;
opts.use_ro_glc      =  false;

% define the thermal and spectral grid
opts.dz_thermal      =  0.04;       % dz for heat transfer          [m]
opts.dz_spectral     =  0.002;      % dz for radiative transfer     [m]
opts.z0_thermal      =  12; %20;         % thickness for heat transfer   [m]
opts.z0_spectral     =  4; %8;          % thickness for rad transfer    [m]
opts.f_ice_min       =  0.01;

opts.ro_snow_i       =  900.0;
opts.i_grainradius   =  25;         % index 25 = 2.0 mm             [#]
opts.nwavl           =  118;
opts.nradii          =  35;

% solver options
opts.fzero           =  optimset('Display','off','TolX',1e-6);

% model structure
if strcmp(meltmodel,'icemodel')
   opts.icemodel  = true;
   opts.skinmodel = false;
elseif strcmp(meltmodel,'skinmodel')
   opts.icemodel  = false;
   opts.skinmodel = true;
end

% forcing data
opts.forcingData     =  forcingdata;
opts.forcingUserData =  userdata;
opts.forcingUserVars =  uservars;

% let METINIT know whether userData is requested
if strcmp(opts.forcingUserData,'none')
   opts.useUserData  =  false;
else
   opts.useUserData  =  true;
end


% to use forcing data from nearby met-stations for catchment runs, the met
% file needs the met station name instead of the catchment 'siteName'

metname = sitename;     % default is siteName

% swap out the met station name if requested
if contains(forcingdata,'KANL') && strcmp(sitename,'upperBasin')
   metname = 'KANL';
elseif strcmp(forcingdata,'KANM') && strcmp(sitename,'behar')
   metname = 'KANM';
elseif strcmp(forcingdata,'KANM') && strcmp(sitename,'slv1')
   metname = 'KANM';
elseif strcmp(forcingdata,'KANM') && strcmp(sitename,'slv2')
   metname = 'KANM';
end

% this sets the met filename. note: for kanm/kanl, all values are from
% the weather station. otherwise, the values are from mar
switch opts.dt
   case 900
      opts.metfname = ['met_' metname '_' forcingdata '_' simyear '_15m.mat'];
   case 3600
      opts.metfname = ['met_' metname '_' forcingdata '_' simyear '_1hr.mat'];
end

% these are experimental options that should be left as-is
opts.tlagcolumn = 6*3600/opts.dt;
opts.tlagsurf   = 12;

% build a filename string to save the output data
fsave = [meltmodel '_' sitename '_' simyear '_' forcingdata '_swap_' 
   upper(userdata) '_' uservars{1}];

opts.fsave = fsave;

opts.error = '';

% temporary until i figure out how/if to combine region runs with site
% runs, or allow for multi-year site runs
opts.simyears  = startyear:1:endyear;
opts.numyears  = endyear-startyear+1;

%    if opts.use_init == true
%       opts.f_init = [opts.path.initData 'v10/run1/init.mat'];
%    end

% % it doesn't make sense to add ALL the vars, if we ever swap more than
% one, but if so, need to remove '_' userVars{1} from above fsave
% % concatenate the swapped variable(s)
%    for mm = 1:numel(unique(userVars))
%       fsave    = [fsave '_' userVars{mm}];
%    end
%    opts.fsave  = fsave;


%     % this is just a helper b/c sometimes I forget to set the year
%     switch siteName
%         case 'upperBasin'
%             simYear     =   '2016';
%         case {'slv1','slv2'}
%             simYear     =   '2015';
%         otherwise
%     end
