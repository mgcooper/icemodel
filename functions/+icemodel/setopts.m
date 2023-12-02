function opts = setopts(simmodel, sitename, simyears, forcings, ...
      userdata, uservars, savedata, casename, testname)
   %SETOPTS Set model options
   %
   %
   % See also: icemodel.config

   if nargin < 8 || isempty(casename); casename = ''; end
   if nargin < 9 || isempty(testname); testname = ''; end

   % Set the project-level configuration
   icemodel.config();
   
   %---------------------------- save the standard options that were passed in
   %----------------------------------------------------------------------------
   opts.savedata = savedata;
   opts.simmodel = simmodel;
   opts.sitename = sitename;
   opts.forcings = forcings;
   opts.userdata = userdata;
   opts.uservars = uservars;
   opts.simyears = simyears;
   opts.numyears = numel(simyears);
   opts.casename = casename;
   opts.testname = testname;

   %---------------------------- optional settings / parameters
   %----------------------------------------------------------------------------

   % general model settings
   opts.spinup_loops    =  1;       % number of spin-up loops to initialize
   opts.use_init        =  false;   % use pre-initialized data?
   opts.kabs_user       =  true;    % use user-defined ice absorptivity?
   opts.use_ro_glc      =  false;   % use same density for liquid/solid ice?
   opts.calendar_type   =  'noleap';

   % model parameters
   opts.z_0             =  0.001;   % Surface aero. roughness length    [m]
   opts.z_obs           =  3.0;     % Height of wind and air temp obs   [m]
   opts.ro_snow_i       =  900.0;   % initial ice density               [kg/m3]
   opts.liqresid        =  0.07;    % residual pore water fraction      [-]

   % timestepping / grid thickness
   opts.dt              =  900.0;   % timestep                             [s]
   opts.dz_thermal      =  0.04;    % dz for heat transfer                 [m]
   opts.dz_spectral     =  0.002;   % dz for radiative transfer            [m]
   opts.z0_thermal      =  20;      % domain thickness for heat transfer   [m]
   opts.z0_spectral     =  4;       % domain thickness for rad transfer    [m]
   opts.f_ice_min       =  0.01;

   % solver options
   opts.sebsolver = 1;
   opts.fzero = optimset('Display', 'off', 'TolX', 1e-6);

   % the mie scattering coefficients are defined for 35 grain sizes and 118
   % spectral bands. define those dimensions here, they are used to read in
   % the data array in GETSCATTERCOEFS. also set the grain size index.
   opts.nwavl           =  118;
   opts.nradii          =  35;
   opts.i_grainradius   =  25;      % index 25 = 2.0 mm                 [#]

   % define the timescale beyond which stored meltwater is assumed to runoff
   opts.tlagcolumn = 6 * 3600 / opts.dt; % convert hours to timesteps
   opts.tlagsurf = 6 * 3600 / opts.dt;
   opts.error = '';

   %---------------------------- set the input and output paths
   %----------------------------------------------------------------------------

   % WRITEOUTPUT appends ['ice1_' opts.casename '.mat'] and saves the file in
   % a subfolder of opts.pathoutput for each year e.g. opts.pathoutput/2016. 
   % 
   % For sector runs, opts.casename is set to the grid point ID outside this
   % function in a loop, to get ice1_1.mat, ice1_2.mat, and so forth. Same for
   % metfname. When writing to netcdf is implemented, that issue will go away.

   opts.pathinput = getenv('ICEMODELINPUTPATH');
   opts.pathoutput = fullfile(getenv('ICEMODELOUTPUTPATH'), sitename, simmodel);

   assert(isfolder(opts.pathinput), ...
      'ICEMODELINPUTPATH does not exist, set it using icemodel.config');

   % For test runs, option to create a subfolder in ICEMODELOUTPUTPATH
   if ~strcmp(opts.testname, 'none')
      opts.pathoutput = fullfile(opts.pathoutput, testname);
   end

   % Create the casename. WRITEOUTPUT appends this to the base filenames.
   opts.casename = icemodel.setcase(forcings, userdata, uservars);
      
   % Create folders for each simulation year in pathoutput/ if they don't exist
   if opts.savedata
      
      icemodel.mkfolders(opts);

      % save the model opts
      optsfile = ['opts_' opts.casename '.mat'];
      save(fullfile(opts.pathoutput, 'opts', optsfile), 'opts');
   end

   %---------------------------- set the met forcing file name
   %----------------------------------------------------------------------------

   if strcmp(sitename, 'sector')

      % Could add logic here to deal with sector file names. For now, the
      % metfname must be set outside this function in a loop.
      % for n = 1:numel(runpoints)
      %    opts.metfname = 'met_sector.mat';
      % end

      opts.vars1 = {'Tsfc'};
      opts.vars2 = {'Tice', 'f_ice', 'f_liq', 'df_liq', 'df_drn', 'df_evp'};

   else

      opts.vars1 = ...
         {'Tsfc', 'Qm', 'Qf', 'Qe', 'Qh', 'Qc', 'chi', 'balance', 'dt_sum'};
      opts.vars2 = ...
         {'Tice', 'f_ice', 'f_liq', 'df_liq', 'df_drn', 'df_evp', 'Sc', ...
         'errH', 'errT', 'lcflag'};

      if strcmp(simmodel, 'skinmodel')
         opts.vars2 = {'Tice', 'f_ice', 'f_liq'};
      end

      % Create met file names
      opts = createMetFileNames(opts, sitename, forcings, simyears);

      % append the path to the met file and out file names
      opts.metfname = fullfile(opts.pathinput, 'met', opts.metfname);

      % % This is the previous method.
      % for n = 1:numel(simyears)
      %    simyear = num2str(simyears(n));
      %    opts.metfname{n} = ['met_' metname '_' forcings '_' simyear '_' dtstr];
      %    opts.outfname{n} = [simmodel '_' sitename '_' simyear '_' forcings ...
      %       '_swap_' upper(userdata) '_' uservars];
      % end
   end
end

%% Create met file names
function opts = createMetFileNames(opts, sitename, forcings, simyears)

   % Deal with the case where met-station forcing data (as opposed to gridded
   % climate model forcing data) is requested for a nearby catchment by
   % replacing the catchment name in the metfile with the met station name.
   % For example, if sitename=="behar" and forcingdata=="kanm", this sets the
   % metfile name to met_kanm_kanm_YYYY rather than met_behar_kanm_YYYY, to
   % negate the need to create a second (identical) met_behar_kanm_YYYY file.
   if strcmpi(forcings, 'kanl') ...
         && ismember(sitename, {'ak4','upperbasin'})
      metname = 'kanl';
   elseif strcmpi(forcings, 'kanm') ...
         && ismember(sitename, {'slv1','slv2','behar'})
      metname = 'kanm';
   else
      metname = sitename;
   end

   % Deal with the option to use 15 min vs 1 hrly forcings.
   switch opts.dt
      case 900
         dtstr = '15m.mat';
      case 3600
         dtstr = '1hr.mat';
   end

   % Build full file names.
   for n = 1:numel(simyears)
      simyear = num2str(simyears(n));
      opts.metfname{n} = ['met_' metname '_' forcings '_' simyear '_' dtstr];
   end
end
