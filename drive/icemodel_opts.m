function opts = icemodel_opts(sitename, simmodel, simyears, forcings, ...
   userdata, uservars, savedata, testname, testpoint)
%ICEMODEL_OPTS set model options

if nargin < 8 || isempty(testname); testname = 'none'; end
if nargin < 9 || isempty(testpoint); testpoint = 'none'; end

%---------------------------- set the standard options that were passed in
%------------------------------------------------------------------------------
opts.savedata = savedata;
opts.simmodel = simmodel;
opts.sitename = sitename;
opts.forcings = forcings;
opts.userdata = userdata;
opts.uservars = uservars;
opts.simyears = simyears;
opts.numyears = numel(simyears);
opts.testname = testname;
opts.testpoint = testpoint;

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
opts.dt              =  900.0;   % timestep                             [s]
opts.dz_thermal      =  0.04;    % dz for heat transfer                 [m]
opts.dz_spectral     =  0.002;   % dz for radiative transfer            [m]
opts.z0_thermal      =  20;      % domain thickness for heat transfer   [m]
opts.z0_spectral     =  4;       % domain thickness for rad transfer    [m]
opts.f_ice_min       =  0.01;

% solver options
opts.fzero = optimset('Display','off','TolX',1e-6);

% the mie scattering coefficients are defined for 35 grain sizes and 118
% spectral bands. define those dimensions here, they are used to read in
% the data array in GETSCATTERCOEFS. also set the grain size index.
opts.nwavl           =  118;
opts.nradii          =  35;
opts.i_grainradius   =  25;      % index 25 = 2.0 mm                 [#]

%---------------------------- set the met forcing file name
%------------------------------------------------------------------------------

% Deal with the case where met-station forcing data (as opposed to gridded
% climate model forcing data) is requested for a nearby catchment by replacing
% the catchment name in the metfile with the met station name. For example, if
% sitename=="behar" and forcingdata=="kanm", this sets the metfile name to
% met_kanm_kanm_YYYY rather than met_behar_kanm_YYYY, to negate the need to
% create a second (identical) met_behar_kanm_YYYY file.

if strcmpi('kanl', forcings) && any(ismember(sitename,{'ak4','upperbasin'}))
   metname = 'kanl';
elseif strcmpi('kanm', forcings) && any(ismember(sitename,{'slv1','slv2','behar'}))
   metname = 'kanm';
else
   metname = sitename;
end

% Sets the met filenames and output filenames. Note: for kanm/kanl, all values
% are from the weather station. otherwise, the values are from MAR.
switch opts.dt
   case 900
      dtstr = '15m.mat';
   case 3600
      dtstr = '1hr.mat';
end
for n = 1:numel(simyears)
   simyear = num2str(simyears(n));
   opts.metfname{n} = ['met_' metname '_' forcings '_' simyear '_' dtstr];
   opts.outfname{n} = [simmodel '_' sitename '_' simyear '_' forcings ...
      '_swap_' upper(userdata) '_' uservars];
end
% programming note: if swapping multiple vars is supported, will need to extend
% the loop above and add a method to METINIT and/or icemodel

%---------------------------- configure the input and output data paths
%------------------------------------------------------------------------------

% set the input data path
% -----------------------

% give precedence to ICEMODELINPUTPATH before searching elsewhere
if ~isempty(getenv('ICEMODELINPUTPATH'))

   % the ICEMODELIINPUTPATH environment variable is defined
   opts.pathinput = getenv('ICEMODELINPUTPATH');

elseif isfile(opts.metfname{1},'file')

   % the metfile exists somewhere else on the path
   opts.pathinput = strrep(fileparts(which(opts.metfname{1})),'met','');

elseif isfolder(fullfile(pwd, 'input', 'met'))

   % assume the metfile exists in the 'input/met/' directory
   opts.pathinput = fullfile(pwd,'input');

else
   % the metfile does not exist on the path
   % rather than error, use isfield(opts,'metfname') to halt in run script
   opts.msg = 'met file not found';
   return
end

% set the output data path
% ------------------------

if ~isempty(getenv('ICEMODELOUTPUTPATH'))

   % the ICEMODELOUTPUTPATH environment variable is defined, use it
   opts.pathoutput = getenv('ICEMODELOUTPUTPATH');

elseif isfolder(fullfile(pwd, 'output'))

   % an 'output/' directory is present, save the data there
   opts.pathoutput = fullfile(pwd, 'output');

else
   % make a temporary output path
   opts.pathoutput = fullfile(pwd,'icemodel_output_tmp');
   warning(['ICEMODELOUTPUTPATH not found, output path set to ' opts.pathoutput]);
end

% append the path to the met file names and set the userdata path
opts.metfname = fullfile(opts.pathinput, 'met', opts.metfname);
opts.outfname = fullfile(opts.pathoutput, opts.outfname);
opts.userpath = fullfile(opts.pathinput, 'userdata');

% these are experimental options that should be left as-is
% define the timescale beyond which stored meltwater is assumed to runoff
opts.tlagcolumn = 6*3600/opts.dt; % convert hours to timesteps
opts.tlagsurf = 6*3600/opts.dt;
opts.error = '';

