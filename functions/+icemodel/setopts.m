function opts = setopts(simmodel, sitename, simyears, forcings, ...
      userdata, uservars, savedata, casename, testname)
   %SETOPTS Set model options
   %
   %
   % See also: icemodel.config

   if nargin < 8 || isempty(casename); casename = sitename; end
   if nargin < 9 || isempty(testname); testname = 'none'; end

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

   % timestepping / grid thickness
   opts.dt              =  900.0;   % timestep                             [s]
   opts.dz_thermal      =  0.04;    % dz for heat transfer                 [m]
   opts.dz_spectral     =  0.002;   % dz for radiative transfer            [m]
   opts.z0_thermal      =  20;      % domain thickness for heat transfer   [m]
   opts.z0_spectral     =  4;       % domain thickness for rad transfer    [m]
   opts.f_ice_min       =  0.01;

   % solver options
   opts.fzero = optimset('Display', 'off', 'TolX', 1e-6);

   % the mie scattering coefficients are defined for 35 grain sizes and 118
   % spectral bands. define those dimensions here, they are used to read in
   % the data array in GETSCATTERCOEFS. also set the grain size index.
   opts.nwavl           =  118;
   opts.nradii          =  35;
   opts.i_grainradius   =  25;      % index 25 = 2.0 mm                 [#]

   % define the timescale beyond which stored meltwater is assumed to runoff
   opts.tlagcolumn = 6*3600/opts.dt; % convert hours to timesteps
   opts.tlagsurf = 6*3600/opts.dt;
   opts.error = '';

   %---------------------------- set the input and output paths
   %----------------------------------------------------------------------------
   opts.pathinput = getenv('ICEMODELINPUTPATH');
   opts.pathoutput = fullfile(getenv('ICEMODELOUTPUTPATH'), ...
      sitename, simmodel, userdata);

   assert(isfolder(opts.pathinput), ...
      'ICEMODELINPUTPATH does not exist, set it using icemodel.config');

   % For test runs, option to create a subfolder in ICEMODELOUTPUTPATH
   if ~strcmp(opts.testname, 'none')
      opts.pathoutput = fullfile(opts.pathoutput, testname);
   end

   % build the output folders if they don't already exist. This creates one
   % folder for each simulation year in pathoutput/
   if opts.savedata
      if ~isfolder(opts.pathoutput)
         warning( ...
            'output folders do not exist, creating them in %s', opts.pathoutput)
         icemodel.mkfolders(opts);
      end

      % Between here was in region not point
      % make a folder to save the model options
      if ~isfolder(fullfile(opts.pathoutput, 'opts'))
         mkdir(fullfile(opts.pathoutput, 'opts'));
      end

      % save the model opts
      optsfile = ['opts_' simmodel '_' userdata '.mat'];
      save(fullfile(opts.pathoutput, 'opts', optsfile), 'opts');
      % Between here was in region not point
   end

   %---------------------------- set the met forcing file name
   %----------------------------------------------------------------------------

   % to reconcile sector and site-scale runs, could try to put the common logic
   % for met file and output file names in the first part of the file name and
   % use empty chars for the non-common or an if-else. For now, set the metfname
   % and outfname in the run script for sector scale runs otherwise

   % opts.metfname = fullfile( ...
   %    opts.pathinput, 'met', ['met_' metname '_<other parts>.mat']);
   % opts.outfname = fullfile( ...
   %    opts.pathoutput, [simmodel '_' sitename '_<other parts>']);

   if strcmp(sitename, 'sector')

      % Could add logic here to deal with sector file names
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
         'errH', 'errT'};

      % Deal with the case where met-station forcing data (as opposed to gridded
      % climate model forcing data) is requested for a nearby catchment by
      % replacing the catchment name in the metfile with the met station name.
      % For example, if sitename=="behar" and forcingdata=="kanm", this sets the
      % metfile name to met_kanm_kanm_YYYY rather than met_behar_kanm_YYYY, to
      % negate the need to create a second (identical) met_behar_kanm_YYYY file.

      if strcmpi(forcings, 'kanl') && ...
            ismember(sitename, {'ak4','upperbasin'})
         metname = 'kanl';
      elseif strcmpi(forcings, 'kanm') && ...
            ismember(sitename, {'slv1','slv2','behar'})
         metname = 'kanm';
      else
         metname = sitename;
      end

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
      % append the path to the met file and out file names
      opts.metfname = fullfile(opts.pathinput, 'met', opts.metfname);
      opts.outfname = fullfile(opts.pathoutput, opts.outfname);
   end
end
