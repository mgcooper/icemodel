function opts = setopts(smbmodel, sitename, simyears, forcings, ...
      userdata, uservars, saveflag, testname, backupflag)
   %SETOPTS Set model options
   %
   %
   % See also: icemodel.config

   if nargin < 5  || isempty(userdata); userdata = 'none'; end
   if nargin < 6  || isempty(uservars); uservars = 'albedo'; end
   if nargin < 7  || isempty(saveflag); saveflag = false; end
   if nargin < 8  || isempty(testname); testname = ''; end
   if nargin < 9  || isempty(backupflag); backupflag = true; end

   [smbmodel, sitename, forcings, userdata, uservars, testname] ...
      = convertStringsToChars(...
      smbmodel, sitename, forcings, userdata, uservars, testname);

   %------------------------- save the standard options passed in
   %--------------------------------------------------------------
   opts.saveflag = saveflag;
   opts.smbmodel = smbmodel;
   opts.sitename = sitename;
   opts.forcings = forcings;
   opts.userdata = userdata;
   opts.uservars = uservars;
   opts.simyears = simyears;
   opts.numyears = numel(simyears);
   opts.testname = testname;
   opts.backupflag = backupflag;

   %------------------------- optional settings / parameters
   %---------------------------------------------------------

   % general model settings
   opts.spinup_loops    =  1;       % number of spin-up loops to initialize
   opts.use_init        =  false;   % use pre-initialized data?
   opts.kabs_user       =  true;    % use user-defined ice absorptivity?
   opts.use_ro_glc      =  false;   % use same density for liquid/solid ice?
   opts.calendar_type   =  'noleap';

   % model parameters
   opts.z_0             =  0.001;   % Surface aero. roughness length    [m]
   opts.ro_ice_init     =  900.0;   % initial ice density               [kg/m3]
   opts.T_ice_init      = -8.0;     % initial ice temperature           [C]
   opts.f_liq_resid     =  0.02;    % residual pore water fraction      [-]

   % solver options and timestepping / grid thickness
   if strcmp(smbmodel, 'icemodel')

      opts.seb_solver      = 1;     % recommended: 1 (1=analytic, 2=numeric)
      opts.bc_type         = 1;     % recommended: 1 (1=dirichlet, 2=robin)
      opts.conduct_type    = 1;     % recommended: 1 (Patankar practice "B")
      opts.maxiter         = 50;    % 1d nonlinear heat transfer max iterations

      opts.dt              = 3600;  % timestep (3600 or 900)               [s]
      opts.dz_thermal      = 0.04;  % dz for thermal heat transfer         [m]
      opts.dz_spectral     = 0.002; % dz for radiative heat transfer       [m]
      opts.z0_thermal      = 20;    % domain thickness for heat transfer   [m]
      opts.z0_spectral     = 8;     % domain thickness for rad transfer    [m]
      opts.f_ice_min       = 0.1;   % layer combination threshold (ice fraction)

   elseif strcmp(smbmodel, 'skinmodel')

      opts.seb_solver      = 1;     % recommended: 1 (1=analytic, 2=numeric)
      opts.bc_type         = 1;     % recommended: 1 (1=dirichlet, 2=robin)
      opts.conduct_type    = 1;     % recommended: 1 (Patankar practice "B")
      opts.maxiter         = 50;    % 1d nonlinear heat transfer max iterations

      opts.dt              = 3600;  % timestep (3600 or 900)               [s]
      opts.dz_thermal      = 0.04;  % dz for thermal heat transfer         [m]
      opts.dz_spectral     = 0.002; % dz for radiative heat transfer       [m]
      opts.z0_thermal      = 12;    % domain thickness for heat transfer   [m]
      opts.z0_spectral     = 4;     % domain thickness for rad transfer    [m]
      opts.f_ice_min       = 0.1;   % layer combination threshold (ice fraction)
   else
      error('unrecognized surface mass balance model name SMBMODEL')
   end

   % The mie scattering coefficients are defined for 35 grain sizes and 118
   % spectral bands. Define those dimensions here, they are used to read in
   % the data array in GETSCATTERCOEFS. Also set the grain size index.
   opts.nwavl           = 118;
   opts.nradii          = 35;
   opts.i_grainradius   = 25;          % index 25 = 2.0 mm              [#]

   % Define the height of the air temperature and wind observations     [m]
   if strcmp(forcings, 'mar')

      opts.z_tair = 2.0;
      opts.z_wind = 10.0;

   elseif any(strcmp(forcings, {'kanl', 'kanm'}))

      opts.z_tair =  2.0; % 2.5
      opts.z_wind =  2.0; % 3.0

   else
      % USER DEFINED
      opts.z_tair =  2.0;
      opts.z_wind =  3.0;
   end

   %------------------------- End of user-defined model options
   %--------------------------------------------------------------

   opts = setRunConfiguration(opts); % Run the automated configuration

   function opts = setRunConfiguration(opts)

      %------------------------- set the input and output paths
      %---------------------------------------------------------

      % WRITEOUTPUT appends ['ice1_' opts.casename '.mat'] and saves the file in
      % a subfolder of opts.pathoutput for each year e.g. opts.pathoutput/2016.
      %
      % For gridded sector runs, opts.casename is set to the grid point ID
      % outside this function in a loop, to get ice1_1.mat, ice1_2.mat, and so
      % forth. Same for metfname. Writing to netcdf should eliminate this.

      opts.pathinput = getenv('ICEMODEL_INPUT_PATH');
      opts.pathoutput = fullfile( ...
         getenv('ICEMODEL_OUTPUT_PATH'), sitename, smbmodel);

      if strcmp(sitename, 'sector')
         if strcmp(userdata, 'none')
            % If userdata is "none", prevent creation of a "none" folder, use
            % the forcings instead e.g. sector/icemodel/mar.
            opts.pathoutput = fullfile(opts.pathoutput, forcings);
         else
            % e.g. sector/icemodel/modis.
            opts.pathoutput = fullfile(opts.pathoutput, userdata);
         end
      end

      assert(isfolder(opts.pathinput), ...
         'ICEMODEL_INPUT_PATH does not exist, set it using icemodel.config');

      % For test runs, option to create a subfolder in ICEMODEL_OUTPUT_PATH
      if ~strcmp(opts.testname, 'none')
         opts.pathoutput = fullfile(opts.pathoutput, testname);
      end

      % Create the casename. WRITEOUTPUT appends this to the base filenames.
      opts.casename = icemodel.setcase(forcings, userdata, uservars);

      % Create folders for each simulation year in output/ if they don't exist.
      if saveflag

         icemodel.mkfolders(opts);

         % save the model opts
         optsfile = fullfile( ...
            opts.pathoutput, 'opts', ['opts_' opts.casename '.mat']);
         backupfile(optsfile, backupflag);
         save(optsfile, 'opts');
      end

      %---------------------------- set the met forcing file name
      %----------------------------------------------------------------------------

      if strcmp(sitename, 'sector')

         % Could add logic here to deal with sector file names. For now, the
         % metfname must be set outside this function in a loop.
         % for n = 1:numel(runpoints)
         %    opts.metfname = 'met_sector.mat';
         % end
      else
         % Create met file names
         opts = createMetFileNames(opts, sitename, forcings, simyears);

         % append the path to the met file and out file names
         opts.metfname = fullfile(opts.pathinput, 'met', opts.metfname);
      end

      %---------------------- set the output file variables
      %-----------------------------------------------------

      % ice1.vars1 = Surface (1-d) data
      % ice2.vars2 = Subsurface (2-d) data

      % Two important programming notes:
      %
      % 1) The order in which the opts.vars1 and opts.vars2 variables are set
      % must match the order in which the data is stored in the cell arrays that
      % are passed to SAVEOUTPUT from icemodel main.
      %
      % 2) If new variables are added, POSTPROC must be reviewed to ensure
      % correct processing is applied, including rounding precision.

      if strcmp(sitename, 'sector')

         if strcmp(smbmodel, 'skinmodel')
            opts.vars1 = {'Tsfc', 'Qm', 'Qe'};
            opts.vars2 = {'Tice'};
         else

            opts.vars1 = {'Tsfc'};
            opts.vars2 = {'Tice', 'f_ice', 'f_liq', 'df_liq', 'df_evp'};
         end

      else

         if strcmp(smbmodel, 'skinmodel')
            % opts.vars1 = {'Tsfc', 'Qm', 'Qe'};
            % opts.vars2 = {'Tice'};
            % % opts.vars2 = {'Tice', 'f_ice', 'f_liq'};

            opts.vars1 = ...
               {'Tsfc', 'Qm', 'Qe', 'Qh', 'Qc', 'chi', 'balance', ...
               'dt_sum', 'Tsfc_converged', 'Tice_converged', 'Tice_numiter'};

            opts.vars2 = ...
               {'Tice', 'f_ice', 'f_liq'};

         else

            opts.vars1 = ...
               {'Tsfc', 'Qm', 'Qe', 'Qh', 'Qc', 'chi', 'balance', ...
               'dt_sum', 'Tsfc_converged', 'Tice_converged', 'Tice_numiter'};

            opts.vars2 = ...
               {'Tice', 'f_ice', 'f_liq', 'df_liq', 'df_evp', 'df_lyr', 'Sc'};
         end
         % % model diagnostics
         % opts.vars3 = ...
         %    {'Tsfc_converged', 'Tsfc_numiter', 'Tice_converged', ...
         %    'Tice_numiter', 'dt_sum', 'Enthalpy_residual', 'Tice_residual', ...
         %    'Layer_combo'};
      end
   end
end
function opts = createMetFileNames(opts, sitename, forcings, simyears)
   %CREATEMETFILENAMES Create icemodel met file names for model opts struct
   %
   %  opts = createMetFileNames(opts, sitename, forcings, simyears)
   %
   % See also: icemodel.setopts

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
   opts.tlag = 6 * 3600 / opts.dt; % hours to timesteps

   % Build full file names.
   for n = 1:numel(simyears)
      simyear = num2str(simyears(n));
      opts.metfname{n} = ['met_' metname '_' forcings '_' simyear '_' dtstr];
   end
end
