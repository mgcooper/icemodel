function opts = setopts(smbmodel, sitename, simyears, forcings, ...
      userdata, uservars, testname, saveflag, backupflag)
   %SETOPTS Set model options
   %
   %  opts = setopts(smbmodel, sitename, simyears, forcings, ...
   %     userdata, uservars, testname, saveflag, backupflag)
   %
   % Inputs
   %
   %  SMBMODEL - "surface mass balance" model option. Options are "icemodel" or
   %  "skinmodel". "icemodel" runs a subsurface energy balance with meltwater
   %  storage and refreezing. "skinmodel" runs a "skin" surface energy balance
   %  with meltwater production but no storage or refreezing.
   %
   %  SITENAME - a string scalar indicating the site name which is used to
   %  locate the forcing data file and create the output file names.
   %
   %  SIMYEARS - a numeric scalar or vector of 4-digit years e.g. 2016, also
   %  used to locate the correct forcing data file.
   %
   %  FORCINGS - a string scalar indicating the forcing data. Options are "mar",
   %  "merra", and "racmo", indicating one of three climate models from which
   %  icemodel forcing files have been created. Users who wish to add new
   %  forcing files can update the icemodel.completions and arguments blocks
   %  used throughout the repository.
   %
   %  USERDATA - (optional) a string scalar indicating the alternative forcing
   %  data to be used in place of the forcing data in the forcing met file.
   %  Options are "mar", "merra", "racmo", "modis", "kanl", and "kanm",
   %  indicating three climate models, modis satellite albedo, and the kanl and
   %  kanm weather stations. The default value is the value of the FORCINGS
   %  parameter.
   %
   %  USERVARS - (optional) a string scalar indicating the variable name in the
   %  userdata forcing file to be "swapped out" with the corresponding column in
   %  the met forcing file. Any variable in the met file which also exists in
   %  the userdata file can be swapped out, but only one variable can be swapped
   %  out at a time.
   %
   %  TESTNAME - (optional) a string scalar indicating a unique run ID. If
   %  supplied, an additional output folder will be nested under the
   %  ICEMODEL_OUTPUT_PATH/SITENAME/SMBMODEL default output folder structure.
   %  This can be used to create an ensemble of model runs e.g.:
   %     ICEMODEL_OUTPUT_PATH/SITENAME/SMBMODEL/run001
   %     ICEMODEL_OUTPUT_PATH/SITENAME/SMBMODEL/run002
   %     ...
   %     ICEMODEL_OUTPUT_PATH/SITENAME/SMBMODEL/run00N
   %
   %  SAVEFLAG - (optional) a logical scalar indicating whether to save the
   %  output data or not.
   %
   %  BACKUPFLAG - (optional) a logical scalar indicating whether to backup the
   %  output files if they already exist.
   %
   % Outputs
   %
   %  OPTS - a struct containing the runtime model options.
   %
   % Boundary-condition / coupling notes (icemodel):
   %
   %  "Partitioned" here means the surface SEB and subsurface enthalpy
   %  equations are solved in separate solver blocks (SEBSOLVE/SFCFLIN and
   %  ICEENBAL), not as one monolithic nonlinear system.
   %
   %  bc_type = 1 (Dirichlet with lagged Ts-T coupling):
   %   - Solve the nonlinear SEB for Ts once per full forcing step (SEBSOLVE).
   %   - During that solve, inner iterations use an analytic Newton-Raphson,
   %     numeric "complex step", or derivative-free (Brent's) method, set by
   %     opts.seb_solver (see SEBSOLVE for details). SEBSOLVE may fall back to a
   %     derivative-free root solve if Newton/complex-step fails.
   %   - Outer fixed-point iterations repeatedly update the conductive closure
   %     using the latest Ts iterate, while the subsurface top-node state
   %     (temperature and conductance) remains lagged from the previous accepted
   %     state.
   %   - Hold the converged Ts fixed as the upper boundary condition of the
   %     subsurface enthalpy solver (ICEENBAL) during substeps.
   %   - Classification: partitioned, lagged, weakly coupled across the
   %     surface/subsurface interface (lagged Ts-T coupling at full-timestep
   %     scale, but no substep Ts-T coupling iterations).
   %
   %  bc_type = 2 (Robin with single coupling sweep):
   %   - Use a linearized SEB boundary condition (SFCFLIN) in the subsurface
   %     enthalpy solver (ICEENBAL).
   %   - After each accepted substep, diagnose Ts from the updated top-node
   %     state and refresh the SEB linearization for the next substep.
   %   - No inner Ts-T convergence loop within a substep (single coupling
   %     sweep per substep).
   %   - Classification: partitioned, lagged, weakly coupled.
   %
   %  bc_type = 3 (Robin with strong Ts-T coupling iterations):
   %   - Within each substep, perform outer fixed-point coupling iterations
   %     between the SEB linearization and subsurface enthalpy solve until Ts
   %     converges (with relaxation).
   %   - Typically slower than bc_type = 2; cost depends on cpltol, maxcpliter,
   %     relaxation, and timestep adaptation.
   %   - Classification: partitioned, strongly coupled at substep scale
   %     (iterative block/Picard coupling, not monolithic Newton).
   %
   % See also: icemodel.config icemodel.run.point

   % Parse inputs
   narginchk(4, 9)

   if nargin < 5  || isempty(userdata); userdata = forcings; end
   if nargin < 6  || isempty(uservars); uservars = 'albedo'; end
   if nargin < 7  || isempty(testname); testname = ''; end
   if nargin < 8  || isempty(saveflag); saveflag = false; end
   if nargin < 9  || isempty(backupflag); backupflag = true; end

   % Legacy option used 'none', set to forcings instead
   if strcmp(userdata, 'none')
      userdata = forcings;
   end
   if strcmp(testname, 'none')
      testname = '';
   end

   % convertStringsToChars in a pre-R2017b compatible way:
   args = {smbmodel, sitename, forcings, userdata, uservars, testname};
   for n = 1:numel(args)
      if isstring(args{n})
         args{n} = char(args{n});
      end
   end
   [smbmodel, sitename, forcings, userdata, uservars, testname] = deal(args{:});

   % If >= R2017b:
   % [smbmodel, sitename, forcings, userdata, uservars, testname] ...
   %    = convertStringsToChars(...
   %    smbmodel, sitename, forcings, userdata, uservars, testname);

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

   % solver, timestepping, and mesh options
   if strcmp(smbmodel, 'icemodel')

      % Solver options. See function doc for info about each bc type.

      opts.bc_type         = 2;     % recommended: 1 (1=dirichlet, 2/3=robin)
      opts.seb_solver      = 1;     % recommended: 1 (1=analytic, 2=numeric)
      opts.conduct_type    = 1;     % recommended: 1 (Patankar practice "B")
      opts.maxiter         = 100;   % inner thermal solver max iterations
      opts.tol             = 1e-2;  % inner thermal solver convergence tolerance [K]
      opts.maxcpliter      = 100;   % outer Ts convergence max iterations
      opts.cpltol          = 1e-2;  % outer Ts convergence tolerance [K]
      opts.sebtol          = 1.0;   % outer SEB convergence tolerance [W m-2]
      opts.omega           = 0.3;   % outer Ts relaxation factor

      opts.dt              = 900;   % timestep (3600 or 900)               [s]
      opts.dz_thermal      = 0.04;  % dz for thermal heat transfer         [m]
      opts.dz_spectral     = 0.002; % dz for radiative heat transfer       [m]
      opts.z0_thermal      = 20;    % domain thickness for heat transfer   [m]
      opts.z0_spectral     = 8;     % domain thickness for rad transfer    [m]
      opts.f_ice_min       = 0.1;   % layer combination threshold (ice fraction)

   elseif strcmp(smbmodel, 'skinmodel')

      opts.bc_type         = 1;     % recommended: 1 (1=dirichlet, 2=robin)
      opts.seb_solver      = 1;     % recommended: 1 (1=analytic, 2=numeric)
      opts.conduct_type    = 1;     % recommended: 1 (Patankar practice "B")
      opts.maxiter         = 100;   % inner thermal solver max iterations
      opts.tol             = 1e-2;  % inner thermal solver convergence tolerance [K]
      opts.maxcpliter      = 100;   % outer Ts convergence max iterations
      opts.cpltol          = 1e-2;  % outer Ts convergence tolerance [K]
      opts.sebtol          = 1.0;   % outer SEB convergence tolerance [W m-2]
      opts.omega           = 1.8;   % outer Ts relaxation factor
      opts.use_aitken      = true;  % use aitken-acceleration or not
      opts.aitken_jumpmax  = 5.0;   % acceleration guess tolerance [K]

      opts.dt              = 900;   % timestep (3600 or 900)               [s]
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

      % Custom option for gridded "sector"-scale runs.
      if strcmp(sitename, 'sector')
         % e.g. sector/icemodel/modis.
         opts.pathoutput = fullfile(opts.pathoutput, userdata);
      end

      assert(exist(opts.pathinput, 'dir') == 7, ...
         'ICEMODEL_INPUT_PATH does not exist, set it using icemodel.config');

      % For test runs, option to create a subfolder in ICEMODEL_OUTPUT_PATH
      opts.pathoutput = fullfile(opts.pathoutput, testname);

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
         && ismember(sitename, {'ak4', 'upperbasin'})
      metname = 'kanl';
   elseif strcmpi(forcings, 'kanm') ...
         && ismember(sitename, {'slv1', 'slv2', 'behar'})
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
