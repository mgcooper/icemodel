function opts = icemodel_opts(sitename,meltmodel,forcingdata,userdata, ...
   uservars,startyear,endyear)

%---------------------------- set the standard options that were passed in
%------------------------------------------------------------------------------
opts.sitename        =  sitename;
opts.meltmodel       =  meltmodel;
opts.forcingData     =  forcingdata;
opts.userdata        =  userdata;
opts.uservars        =  uservars;
opts.simyears        =  startyear:1:endyear;
opts.numyears        =  endyear-startyear+1;

%---------------------------- optional settings / parameters
%------------------------------------------------------------------------------

% general model settings
opts.annual_loops    =  1;       % number of spin-up loops to initialize
opts.use_init        =  false;   % use pre-initialized data?
opts.kabs_user       =  true;    % use user-defined ice absorptivity?
opts.use_ro_glc      =  false;   % use same density for liquid/solid ice?
opts.calendar_type   =  'noleap';

% model parameters
opts.z_0             =  0.001;   % Surface aero. roughness length    [m]
opts.z_obs           =  3.0;     % Height of wind and air temp obs   [m]
opts.ro_snow_i       =  900.0;   % initial ice density               [kg/m3]

% timestepping / grid thickness
opts.dt              =  3600.0;   % timestep                             [s]
opts.dz_thermal      =  0.04;    % dz for heat transfer                 [m]
opts.dz_spectral     =  0.002;   % dz for radiative transfer            [m]
opts.z0_thermal      =  12;      % domain thickness for heat transfer   [m]
opts.z0_spectral     =  4;       % domain thickness for rad transfer    [m]
opts.f_ice_min       =  0.01;

% solver options
opts.fzero           =  optimset('Display','off','TolX',1e-6);

% the mie scattering coefficients are defined for 35 grain sizes and 118
% spectral bands. define those dimensions here, they are used to read in
% the data array in GETSCATTERCOEFS. also set the grain size index.
opts.nwavl           =  118;
opts.nradii          =  35;
opts.i_grainradius   =  25;      % index 25 = 2.0 mm                 [#]

% define the timescale beyond which stored meltwater is assumed to runoff
opts.tlagcolumn      =  6*3600/opts.dt;   % convert hours to timesteps
opts.tlagsurf        =  6*3600/opts.dt;

%---------------------------- set the met forcing file name
%------------------------------------------------------------------------------

% deal with the case where met-station forcing data is requested for a
% nearby catchment, as opposed to climate model forcing data that was
% created for the catchment. for example, if sitename=="behar" and
% forcingdata=="kanm", this loads the met_kanm_kanm_YYYY met file, to negate
% the need to create a second (identical) met_behar_kanm_YYYY file.

if forcingdata == "kanl" && any(ismember(sitename,{'ak4','upperbasin'}))
   metname = 'kanl';
elseif forcingdata == "kanm" && any(ismember(sitename,{'slv1','slv2','behar'}))
   metname = 'kanm';
else
   metname = sitename;
end

% multi-year runs are not currently supported, so simyear equals startyear
simyear = num2str(startyear);
switch opts.dt
   case 900
      metfname = ['met_' metname '_' forcingdata '_' simyear '_15m.mat'];
   case 3600
      metfname = ['met_' metname '_' forcingdata '_' simyear '_1hr.mat'];
end

%---------------------------- configure the input and output data paths
%------------------------------------------------------------------------------

% set the input data path
% -----------------------

% give precedence to ICEMODELINPUTPATH before searching elsewhere
if (  ~isempty(getenv('ICEMODELINPUTPATH'))                             ...
      && exist([getenv('ICEMODELINPUTPATH') 'met/' metfname],'file') == 2 )

   % the ICEMODELIINPUTPATH environment variable is defined
   opts.pathinput = getenv('ICEMODELINPUTPATH');

elseif exist([pwd '/input/met/' metfname],'file') == 2

   % the metfile exists in the 'input/met/' directory
   opts.pathinput = [pwd '/input/'];

elseif exist(metfname,'file') == 2

   % the metfile exists somewhere else on the path
   opts.pathinput = strrep(fileparts(which(metfname)),'met','');
else
   % the metfile does not exist on the path
   % rather than error, use isfield(opts,'metfname') to halt in run script
   warning('met file not found');
   return;
end

% set the path to the met file and the user data
opts.metfname = [opts.pathinput 'met/' metfname];
opts.userpath = [opts.pathinput 'userdata/'];

% set the output data path
% ------------------------

if ~isempty(getenv('ICEMODELOUTPUTPATH'))

   % the ICEMODELOUTPUTPATH environment variable is defined, use it
   opts.pathoutput = getenv('ICEMODELOUTPUTPATH');

elseif exist([pwd '/output/'],'dir') == 7

   % an 'output/' directory is present, save the data there
   opts.pathoutput = [pwd '/output/'];

else
   % make a temporary output path
   opts.pathoutput = [pwd '/icemodel_output_tmp/'];
   warning(['ICEMODELOUTPUTPATH not found, output path set to ' opts.pathoutput]);
end

% build a filename string to save the output data
fsave = [ meltmodel '_' sitename '_' simyear '_' upper(forcingdata) ...
   '_swap_' upper(userdata) '_' uservars];

opts.fsave = [opts.pathoutput fsave];




