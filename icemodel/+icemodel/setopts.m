function opts = setopts(smbmodel, sitename, simyears, forcings, ...
      userdata, uservars, testname, saveflag, backupflag, varargin)
   %SETOPTS Set model options
   %
   %  opts = setopts(smbmodel, sitename, simyears, forcings, ...
   %     userdata, uservars, testname, saveflag, backupflag)
   %  opts = setopts(..., 'dt', 900, 'solver', 3)
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
   %  FORCINGS - a string scalar indicating the forcing data. Standard values
   %  include the climate-model forcings "mar", "merra", and "racmo", plus
   %  supported self-forced station runs such as "kanm" and "kanl". Users who
   %  wish to add new forcing files can update the icemodel.namelists
   %  definitions used throughout the repository.
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
   %  Any field in the returned OPTS struct can be overridden by passing
   %  trailing name/value pairs to setopts, or by calling icemodel.resetopts
   %  on the struct before passing it to icemodel or skinmodel. These are the
   %  only two entry points for user-defined parameter values.
   %
   %  The OPTS struct is the complete and authoritative list of user-settable
   %  model parameters. Any parameter not represented in OPTS is fixed at its
   %  canonical value inside icemodel.parameterLookup and cannot be varied
   %  without editing that file directly.
   %
   %  Fields initialized to [] below are user-settable parameters whose
   %  physics-based defaults live in parameterLookup. Leaving them empty
   %  instructs configureRun to fill them from parameterLookup at run time.
   %  Setting a non-empty value (via name/value override or resetopts) bypasses
   %  the lookup and uses the supplied value instead. This is the intended
   %  mechanism for sensitivity analyses and site-specific calibration.
   %
   % Outputs
   %
   %  OPTS - a struct containing the runtime model options.
   %
   % Boundary-condition / coupling notes (icemodel):
   %
   %  "Partitioned" here means the surface SEB and subsurface enthalpy
   %  equations are solved in separate solver blocks
   %  (solve_surface_energy_balance/surface_flux_linearization and
   %  solve_column_enthalpy), not as one monolithic nonlinear system.
   %
   %  "Coupled" here means the two blocks exchange boundary information across
   %  the surface/subsurface interface. "Weakly coupled" means that exchange is
   %  done with one lagged sweep per substep. "Strongly coupled" means the
   %  surface and subsurface blocks are iterated together within the same
   %  substep until the interface state is mutually consistent.
   %
   %  solver = 0 (Dirichlet single-sweep):
   %   - Same algorithm as solver = 1 but with cpl_maxiter forced to 1,
   %     so only one Ts-T coupling pass is performed per substep. Intended
   %     as a lightweight diagnostic mode that exercises the Dirichlet code
   %     path without the inner convergence loop cost. Equivalent behaviour
   %     to solver = 2 but via the Dirichlet coupler rather than Robin.
   %   - Classification: partitioned, weakly coupled (single Dirichlet sweep
   %     per substep, no interface convergence iterations).
   %
   %  solver = 1 (coupled Dirichlet Ts-T iterations):
   %   - Solve the nonlinear SEB for a trial Ts
   %     (icemodel.surface.solve_surface_energy_balance).
   %   - Use that Ts as the upper boundary condition of the subsurface
   %     enthalpy solver (`icemodel.column.solve_column_enthalpy`), then repeat
   %     the surface/subsurface exchange until Ts and the accepted updated-state
   %     SEB residual are mutually consistent within the same substep.
   %   - During each solve_surface_energy_balance call, inner surface
   %     iterations use an analytic Newton-Raphson, numeric "complex step", or
   %     derivative-free (Brent's) method, set by opts.seb_solver (see
   %     solve_surface_energy_balance for details).
   %   - Classification: partitioned and strongly coupled at substep scale
   %     through iterative block/Picard Dirichlet coupling, not monolithic
   %     Newton over the full surface-subsurface state.
   %
   %  solver = 2 (Robin with single sweep):
   %   - Use a linearized SEB boundary condition
   %     (surface_flux_linearization) in the subsurface enthalpy solver
   %     (`icemodel.column.solve_column_enthalpy`).
   %   - After each accepted substep, diagnose Ts from the updated top-node
   %     state and refresh the SEB linearization for the next substep.
   %   - No inner Ts-T convergence loop within a substep (single coupling
   %     sweep per substep).
   %   - Note: this is the one-iteration special case of the same partitioned
   %     Robin coupler used by solver = 3, implemented by forcing cpl_maxiter=1
   %     so each substep performs exactly one Robin sweep.
   %   - Classification: partitioned, lagged, weakly coupled (one explicit
   %     Robin sweep per substep, no interface convergence iterations).
   %
   %  solver = 3 (Robin with strong Ts-T coupling iterations):
   %   - Within each substep, perform outer fixed-point coupling iterations
   %     between the SEB linearization and subsurface enthalpy solve until Ts
   %     and SEB residual converge (with relaxation/Aitken acceleration).
   %   - Typically slower than solver = 2; cost depends on cpl_Ts_tol,
   %     cpl_maxiter, relaxation, and timestep adaptation.
   %   - Classification: partitioned, strongly coupled at substep scale
   %     (iterative block/Picard coupling, not monolithic Newton).
   %
   % See also: icemodel.config icemodel.run.point icemodel.configureRun
   % icemodel.getopts icemodel.resetopts

   % Parse inputs.
   % Keep explicit positional parsing here instead of an arguments block for
   % pre-R2019b compatibility. Namespace/test helpers may use arguments
   % blocks, but this core runtime entry point should preserve the older
   % calling contract until that compatibility target is intentionally dropped.
   narginchk(4, inf)
   if nargin < 5; userdata = []; end
   if nargin < 6; uservars = []; end
   if nargin < 7; testname = []; end
   if nargin < 8; saveflag = []; end
   if nargin < 9; backupflag = []; end
   [smbmodel, sitename, simyears, forcings, userdata, uservars, ...
      testname, saveflag, backupflag] = parseinputs(nargin, ...
      smbmodel, sitename, simyears, forcings, userdata, uservars, ...
      testname, saveflag, backupflag);

   %------------------------- save the standard options passed in
   %--------------------------------------------------------------
   opts = initopts(smbmodel, sitename, simyears, forcings, ...
      userdata, uservars, testname, saveflag, backupflag);

   %------------------------- optional settings / parameters
   %---------------------------------------------------------

   %%% General model settings
   opts.n_spinup_years  =  0;       % number of leading simulation years used only for spinup
   opts.use_init        =  false;   % reserved for later generic initialization support
   opts.initfile        =  '';      % reserved generic initialization source
   opts.use_restart     =  false;   % load an exact year-boundary restart state?
   opts.restartfile     =  '';      % restart state file used when use_restart=true
   opts.saverestart     =  false;   % save a restart state at each year boundary

   %%% Debug mode — enable via resetopts(opts, 'debug', true)
   opts.debug           =  false;   % enable solver diagnostic dumps
   opts.debug_path      =  '';      % override for debug output folder

   %%% Model parameters and related options
   opts.use_ro_glc      =  false;   % use same density for liquid/solid ice?
   opts.ro_ice_init     =  900.0;   % initial ice density               [kg/m3]
   opts.T_ice_init      = -8.0;     % initial ice temperature           [C]
   %
   % Capillary residual liquid water per pore volume fractions [-].
   opts.f_res_pore_ice  =  0.02;
   opts.f_res_pore_snow =  0.04;
   opts.f_res_pore_firn =  0.02;

   %%% Timestepping / mesh options
   opts.calendar_type   = 'noleap'; %
   opts.dt              = 900;      % timestep: 3600 or (recommended) 900  [s]
   opts.dz_thermal      = 0.04;     % dz for thermal heat transfer         [m]
   opts.dz_spectral     = 0.002;    % dz for radiative heat transfer       [m]
   opts.z0_thermal      = 20;       % domain thickness for heat transfer   [m]
   opts.z0_spectral     = 8;        % domain thickness for rad transfer    [m]
   opts.f_ice_min       = 0.1;      % minimum ice fraction (remeshing threshold)
   opts.mesh_type       = 1;        % recommended: 1 (Patankar practice "B")

   %%% Surface turbulent-heat-flux scheme.
   %
   % Default = 'bulk_richardson'. 'monin_obukhov' is opt-in via resetopts;
   % its runtime guard is enforced in configureRun.
   opts.turbulent_flux_scheme = 'bulk_richardson';
   %
   % Roughness lengths default to parameterLookup values via configureRun if
   % left empty. Override via resetopts for site-specific or sensitivity runs.
   opts.z0_bulk                = [];
   opts.z0_ice                 = [];
   opts.z0_snow_low_density    = [];
   opts.z0_snow_high_density   = [];
   %
   % Observation heights for the turbulent-flux scheme. Known forcing-dependent
   % defaults are set below after the solver block; Override via resetopts for
   % forcing-specific runs that are not resolved below.
   opts.z_tair = [];   % air temperature observation height [m]
   opts.z_wind = [];   % wind speed observation height      [m]
   opts.z_relh = [];   % relative humidity obs height       [m] (= z_tair)
   %
   % Temporary patch to use forcing snow depth to set roughness lengths until
   % the production snow model is implemented.
   opts.use_forcing_snow_depth_for_thf = false;

   %%% Radiative transfer scheme.
   %
   opts.radiative_transfer_scheme = 'two_stream_schlatter_brandt_warren';
   %
   % The mie coefficients are defined for 118 spectral bands and 35 grain sizes.
   % Define these dimensions for icemodel.radiation.get_scattering_coefficients
   % to read in the data array. Also set the grain size index.
   opts.nwavel          = 118;
   opts.nradii          = 35;
   opts.i_grainradius   = 25; % [#] index 25 = 2.0 mm
   %
   % Option to use a lookup-table to compute bulk extinction coefficients.
   opts.lookup_k_bulk   = true;
   %
   % Option to use a user-defined ice absorptivity spectrum. If true, the model
   % uses the spectrum shipped with the repo from Cooper et al., 2021.
   opts.kabs_user       = true;

   %%% Solver options. See function doc for info about each solver mode.
   if strcmp(smbmodel, 'icemodel')

      % Main solver mode (surface-subsurface coupler)
      % 0 = Dirichlet w/ single Ts-T coupling iteration
      % 1 = Dirichlet w/ strong Ts-T coupling iterations
      % 2 = Robin w/ single Ts-T coupling iteration
      % 3 = Robin w/ strong Ts-T coupling iterations
      opts.solver          = 3;     % recommended: 3

      % Surface (SEB) solver mode. Only relevant for Dirichlet boundary
      % condition (solver = 1). For 'monin_obukhov' thf scheme, seb_solver = 2
      % (numerical Jacobian) is enforced via configureRun.
      opts.seb_solver      = 1;     % recommended: 1 (1=analytic, 2=numeric)

      % Subsurface (column) solver options
      opts.maxiter         = 100;   % thermal solver max iterations
      opts.tol             = 1e-2;  % thermal solver convergence tolerance [K]
      opts.alpha           = 1.0;   % thermal solver relaxation factor (rec: 1.0)
      opts.use_aitken      = false; % thermal solver aitken-acceleration flag
      opts.jumpmax         = 5.0;   % thermal solver acceleration guess tolerance [K]

      % Surface-subsurface coupler options
      opts.cpl_maxiter     = 100;   % coupler Ts convergence max iterations
      opts.cpl_Ts_tol      = 1e-2;  % coupler Ts convergence tolerance [K]
      opts.cpl_seb_tol     = 1.0;   % coupler SEB convergence tolerance [W m-2]
      opts.cpl_alpha       = 1.0;   % coupler Ts relaxation factor (rec: 1.0)
      opts.cpl_aitken      = true;  % coupler Ts aitken-acceleration flag
      opts.cpl_jumpmax     = 5.0;   % coupler Ts acceleration guess tolerance [K]

   elseif strcmp(smbmodel, 'skinmodel')

      % Main solver mode (surface-subsurface coupler)
      opts.solver          = 1;     % required: 1 (2/3 not implemented)

      % Surface (SEB) solver mode
      opts.seb_solver      = 1;     % recommended: 1 (1=analytic, 2=numeric)

      % Subsurface (column) solver options
      opts.maxiter         = 100;   % thermal solver max iterations
      opts.tol             = 1e-2;  % thermal solver convergence tolerance [K]
      opts.alpha           = 1.8;   % thermal solver relaxation factor (rec: 1.8)
      opts.use_aitken      = true;  % thermal solver aitken-acceleration flag (rec: true)
      opts.jumpmax         = 5.0;   % thermal solver acceleration guess tolerance [K]

      % Surface-subsurface coupler options
      opts.cpl_maxiter     = 100;   % coupler Ts convergence max iterations
      opts.cpl_Ts_tol      = 1e-2;  % coupler Ts convergence tolerance [K]
      opts.cpl_seb_tol     = 1.0;   % coupler SEB convergence tolerance [W m-2]
      opts.cpl_alpha       = 1.8;   % coupler Ts relaxation factor (rec: 1.8)
      opts.cpl_aitken      = true;  % coupler Ts aitken-acceleration flag (rec: true)
      opts.cpl_jumpmax     = 5.0;   % coupler Ts acceleration guess tolerance [K]
   else
      error('unrecognized surface mass balance model name SMBMODEL')
   end

   % Set cpl_maxiter=1 for single-sweep Ts-T coupling.
   if ismember(opts.solver, [0 2])
      opts.cpl_maxiter = 1;
   end

   % Forcing-dependent observation heights for the turbulent-flux scheme.
   switch forcings
      case 'mar'
         opts.z_tair = 2.0;
         opts.z_wind = 10.0;

      case 'kanl'
         % Check GEUS/PROMICE metadata, 2.5 and 3.0 may be appropriate.
         opts.z_tair = 2.0;
         opts.z_wind = 2.0;

      case 'kanm'
         % Check GEUS/PROMICE metadata, 2.5 and 3.0 may be appropriate.
         opts.z_tair = 2.0;
         opts.z_wind = 2.0;

      case 'merra'
         opts.z_tair = 2.0;
         opts.z_wind = 2.0;

      otherwise
         % Assume standard heights for other forcings.
         opts.z_tair = 2.0;
         opts.z_wind = 2.0;
   end

   % Assume z_relh = z_tair. This option is not currently used in production.
   opts.z_relh = opts.z_tair;

   % Lag time used by icemodel.column.diagnose_column_runoff, converted from
   % hours to timesteps.
   opts.tlag = 6 * 3600 / opts.dt;

   % Output profile. "minimal" is the lean profile used when wrappers request
   % reduced grid-style output; "standard" is the full point-run profile; and
   % "diagnostic" extends the point-run profile with solver/THF debugging
   % scalars.
   opts.output_profile = 'standard';

   %------------------------- End of user-defined model options
   %--------------------------------------------------------------

   % Initialize run-time contracts. If external callers modify opts via
   % resetopts, contracts are re-enforced by configureRun at model run-time.
   opts = icemodel.resetopts(opts, varargin{:});
   opts = icemodel.configureRun(opts);
end

%%
function [smbmodel, sitename, simyears, forcings, userdata, uservars, ...
      testname, saveflag, backupflag] = parseinputs(nin, smbmodel, ...
      sitename, simyears, forcings, userdata, uservars, testname, ...
      saveflag, backupflag)

   if nin < 5 || isempty(userdata) || isblanktext(userdata)
      userdata = forcings;
   end
   if nin < 6 || isempty(uservars) || isblanktext(uservars)
      uservars = 'albedo';
   end
   if nin < 7 || isempty(testname) || isblanktext(testname)
      testname = '';
   end
   if nin < 8 || isempty(saveflag)
      saveflag = false;
   end
   if nin < 9 || isempty(backupflag)
      backupflag = true;
   end

   if strcmpi(userdata, 'none')
      userdata = forcings;
   end
   if strcmpi(testname, 'none')
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
end

%%
function opts = initopts(smbmodel, sitename, simyears, forcings, ...
      userdata, uservars, testname, saveflag, backupflag)

   opts.saveflag = saveflag;
   opts.smbmodel = smbmodel;
   opts.sitename = sitename;
   opts.forcings = forcings;
   opts.userdata = userdata;
   opts.uservars = uservars;
   opts.simyears = simyears;
   opts.numyears = numel(simyears);
   opts.output_years = [];
   opts.testname = testname;
   opts.saveopts = saveflag;
   opts.backupflag = backupflag;

   % Set defaults for values set in icemodel.configureRun
   opts.pathdata = [];
   opts.pathinput = [];
   opts.patheval = [];
   opts.pathuserdata = [];
   opts.pathoutput = [];
   opts.pathrestart = [];
   opts.casename = [];
   opts.metfname = {};
   opts.vars1 = {};
   opts.vars2 = {};
end
