function opts = setopts(smbmodel, sitename, simyears, forcings, ...
      userdata, uservars, testname, saveflag, backupflag, varargin)
   %SETOPTS Set model options
   %
   %  opts = setopts(smbmodel, sitename, simyears, forcings, ...
   %     userdata, uservars, testname, saveflag, backupflag)
   %  opts = setopts(..., 'dt', 900, 'bc_type', 3)
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
   % Optional name/value overrides
   %
   %  Any field already present in the returned OPTS struct can be overridden
   %  by passing trailing name/value pairs. This is intended for solver,
   %  timestep, and test configuration overrides after the standard run inputs
   %  have been defined by the positional arguments.
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
   %   - Solve the nonlinear SEB for Ts once per substep (SEBSOLVE).
   %   - During that solve, inner iterations use an analytic Newton-Raphson,
   %     numeric "complex step", or derivative-free (Brent's) method, set by
   %     opts.seb_solver (see SEBSOLVE for details).
   %   - Outer fixed-point iterations repeatedly update the conductive closure
   %     using the latest Ts iterate, while the subsurface top-node state
   %     (temperature and conductance) remains lagged from the previous accepted
   %     substep state.
   %   - Use the converged Ts as the upper boundary condition of the subsurface
   %     enthalpy solver (ICEENBAL) for that substep.
   %   - Classification: partitioned, lagged, weakly coupled across the
   %     surface/subsurface interface (one coupling sweep per substep, with
   %     outer iterations that converge Ts against a lagged top-node conductive
   %     closure. No strong Ts–T block-coupling around SEBSOLVE+ICEENBAL within
   %     the substep).
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
   %     between the SEB linearization and subsurface enthalpy solve until
   %     Ts and SEB residual converge (with relaxation/Aitken).
   %   - Typically slower than bc_type = 2; cost depends on cpl_Ts_tol, cpl_maxiter,
   %     relaxation, and timestep adaptation.
   %   - Classification: partitioned, strongly coupled at substep scale
   %     (iterative block/Picard coupling, not monolithic Newton).
   %
   % See also: icemodel.config icemodel.run.point icemodel.configureRun
   % icemodel.getopts icemodel.resetopts

   % Parse inputs
   narginchk(4, inf)

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
   opts.saveopts = saveflag;
   opts.backupflag = backupflag;

   % Set defaults for values set in icemodel.configureRun
   opts.pathinput = [];
   opts.pathoutput = [];
   opts.casename = [];
   opts.metfname = {};
   opts.vars1 = {};
   opts.vars2 = {};

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

      opts.maxiter         = 100;   % thermal solver max iterations
      opts.tol             = 1e-2;  % thermal solver convergence tolerance [K]
      opts.alpha           = 1.0;   % thermal solver relaxation factor (rec: 1.0)
      opts.use_aitken      = false; % thermal solver aitken-acceleration flag
      opts.jumpmax         = 5.0;   % thermal solver acceleration guess tolerance [K]

      opts.cpl_maxiter     = 100;   % coupler Ts convergence max iterations
      opts.cpl_Ts_tol      = 1e-2;  % coupler Ts convergence tolerance [K]
      opts.cpl_seb_tol     = 1.0;   % coupler SEB convergence tolerance [W m-2]
      opts.cpl_alpha       = 1.0;   % coupler Ts relaxation factor (rec: 1.0)
      opts.cpl_aitken      = true;  % coupler Ts aitken-acceleration flag
      opts.cpl_jumpmax     = 5.0;   % coupler Ts acceleration guess tolerance [K]

      % Timestepping / mesh options
      opts.dt              = 900;   % timestep (3600 or (recommended) 900) [s]
      opts.dz_thermal      = 0.04;  % dz for thermal heat transfer         [m]
      opts.dz_spectral     = 0.002; % dz for radiative heat transfer       [m]
      opts.z0_thermal      = 20;    % domain thickness for heat transfer   [m]
      opts.z0_spectral     = 8;     % domain thickness for rad transfer    [m]
      opts.f_ice_min       = 0.1;   % layer combination threshold (ice fraction)

   elseif strcmp(smbmodel, 'skinmodel')

      % Solver options. See function doc for info about each bc type.
      opts.bc_type         = 1;     % recommended: 1 (1=dirichlet, 2=robin)
      opts.seb_solver      = 1;     % recommended: 1 (1=analytic, 2=numeric)
      opts.conduct_type    = 1;     % recommended: 1 (Patankar practice "B")

      opts.maxiter         = 100;   % thermal solver max iterations
      opts.tol             = 1e-2;  % thermal solver convergence tolerance [K]
      opts.alpha           = 1.8;   % thermal solver relaxation factor (rec: 1.8)
      opts.use_aitken      = true;  % thermal solver aitken-acceleration flag (rec: true)
      opts.jumpmax         = 5.0;   % thermal solver acceleration guess tolerance [K]

      opts.cpl_maxiter     = 100;   % coupler Ts convergence max iterations
      opts.cpl_Ts_tol      = 1e-2;  % coupler Ts convergence tolerance [K]
      opts.cpl_seb_tol     = 1.0;   % coupler SEB convergence tolerance [W m-2]
      opts.cpl_alpha       = 1.8;   % coupler Ts relaxation factor (rec: 1.8)
      opts.cpl_aitken      = true;  % coupler Ts aitken-acceleration flag (rec: true)
      opts.cpl_jumpmax     = 5.0;   % coupler Ts acceleration guess tolerance [K]

      % Timestepping / mesh options
      opts.dt              = 900;   % timestep (3600 or (recommended) 900) [s]
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

   % Lag time used by ICERUNOFF, converted from hours to timesteps.
   opts.tlag = 6 * 3600 / opts.dt;

   % Output profile. "minimal" is the lean profile historically tied to
   % gridded sector-scale runs; "standard" is the full point-run profile.
   if strcmp(sitename, 'sector')
      opts.output_profile = 'minimal';
   else
      opts.output_profile = 'standard';
   end

   %------------------------- End of user-defined model options
   %--------------------------------------------------------------

   opts = icemodel.resetopts(opts, varargin{:});
   opts = icemodel.configureRun(opts);
end
