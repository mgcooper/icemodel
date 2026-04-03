function opts = configureRun(opts)
   %CONFIGURERUN Fill derived run settings from OPTS.
   %
   %  opts = icemodel.configureRun(opts)
   %
   % This function is idempotent. It fills required derived fields when they
   % are empty, while preserving caller-supplied values for fields such as
   % PATHINPUT, PATHOUTPUT, CASENAME, METFNAME, VARS1, and VARS2.

   opts = validateTurbulentFluxOptions(opts);

   if ~isfield(opts, 'pathdata') || isempty(opts.pathdata)
      opts.pathdata = icemodel.getpath('data');
   end

   % OUTPUT_YEARS are the post-spinup years retained in saved/postprocessed
   % output. Keep them derived from SIMYEARS and N_SPINUP_YEARS.
   opts.output_years = icemodel.outputYears(opts);

   if ~isfield(opts, 'pathinput') || isempty(opts.pathinput)
      opts.pathinput = icemodel.getpath('input');
   end
   assert(exist(opts.pathinput, 'dir') == 7, ...
      'ICEMODEL_INPUT_PATH does not exist, set it using icemodel.config');

   if ~isfield(opts, 'patheval') || isempty(opts.patheval)
      opts.patheval = icemodel.getpath('eval');
   end

   if ~isfield(opts, 'pathuserdata') || isempty(opts.pathuserdata)
      opts.pathuserdata = icemodel.getpath('userdata');
   end

   % WRITEOUTPUT appends ['ice1_' opts.casename '.mat'] and saves the file in
   % a subfolder of opts.pathoutput for each year e.g. opts.pathoutput/2016.
   %
   % For grid runs, opts.casename and opts.metfname are set outside this
   % function in the wrapper that loops over grid-cell IDs. Core icemodel only
   % provides the canonical base output folder; any legacy extra subfoldering
   % for gridded workflows belongs in the wrapper, not here.
   if ~isfield(opts, 'pathoutput') || isempty(opts.pathoutput)
      opts.pathoutput = icemodel.getpath('output', ...
         opts.sitename, opts.smbmodel, '', [], opts.testname);
   end

   if ~isfield(opts, 'pathrestart') || isempty(opts.pathrestart)
      opts.pathrestart = icemodel.getpath('restart', ...
         opts.sitename, opts.smbmodel, opts.userdata, [], opts.testname);
   end

   % Create the casename. WRITEOUTPUT appends this to the base filenames.
   % For grid runs, the wrapper overwrites this with the grid-cell ID.
   if ~isfield(opts, 'casename') || isempty(opts.casename)
      opts.casename = icemodel.setcase(opts.forcings, opts.userdata, opts.uservars);
   end

   if ~isfield(opts, 'output_profile') || isempty(opts.output_profile)
      opts.output_profile = 'standard';
   end

   if ~isfield(opts, 'metfname') || isempty(opts.metfname)
      opts.metfname = fullfile(opts.pathinput, 'met', ...
         icemodel.createMetFileNames(opts));
   elseif ischar(opts.metfname) || isstring(opts.metfname)
      opts.metfname = cellstr(opts.metfname);
   end

   if ~isfield(opts, 'vars1') || isempty(opts.vars1) ...
         || ~isfield(opts, 'vars2') || isempty(opts.vars2)
      [vars1, vars2] = defaultOutputVariables(opts);
      if ~isfield(opts, 'vars1') || isempty(opts.vars1)
         opts.vars1 = vars1;
      end
      if ~isfield(opts, 'vars2') || isempty(opts.vars2)
         opts.vars2 = vars2;
      end
   end

   % Configure debug diagnostic output when debug mode is enabled via
   % resetopts(opts, 'debug', true). This derives a default debug output
   % folder under ICEMODEL_OUTPUT_PATH/debug/... (paralleling the standard
   % output path structure) and installs the per-kernel environment
   % variables that the solver dump functions already check.
   if isfield(opts, 'debug') && opts.debug
      opts = configureDebugPaths(opts);
   end
end

function opts = validateTurbulentFluxOptions(opts)
   %VALIDATETURBULENTFLUXOPTIONS Validate the turbulent-flux option surface.

   [z0_bulk_default, z0_ice_default, z0_snow_low_default, ...
      z0_snow_high_default] = icemodel.parameterLookup( ...
      'thf_z0_bulk', 'thf_z0_ice', 'thf_z0_snow_low_density', ...
      'thf_z0_snow_high_density');

   if ~isfield(opts, 'turbulent_flux_scheme') ...
         || isempty(opts.turbulent_flux_scheme)
      opts.turbulent_flux_scheme = 'bulk_richardson';
   end

   scheme = lower(char(opts.turbulent_flux_scheme));
   if ~ismember(scheme, {'bulk_richardson', 'bulk_mo'})
      error('icemodel:configureRun:unknownTurbulentFluxScheme', ...
         'Unrecognized turbulent flux scheme: %s', opts.turbulent_flux_scheme);
   end
   opts.turbulent_flux_scheme = scheme;

   if ~isfield(opts, 'z_relh') || isempty(opts.z_relh)
      opts.z_relh = opts.z_tair;
   end

   if ~isfield(opts, 'z0_bulk') || isempty(opts.z0_bulk)
      opts.z0_bulk = z0_bulk_default;
   end
   if ~isfield(opts, 'z0_ice') || isempty(opts.z0_ice)
      opts.z0_ice = z0_ice_default;
   end
   if ~isfield(opts, 'z0_snow_low_density') || isempty(opts.z0_snow_low_density)
      opts.z0_snow_low_density = z0_snow_low_default;
   end
   if ~isfield(opts, 'z0_snow_high_density') || isempty(opts.z0_snow_high_density)
      opts.z0_snow_high_density = z0_snow_high_default;
   end
   if ~isfield(opts, 'use_forcing_snow_depth_for_thf') ...
         || isempty(opts.use_forcing_snow_depth_for_thf)
      opts.use_forcing_snow_depth_for_thf = false;
   end
   opts.use_forcing_snow_depth_for_thf = logical( ...
      opts.use_forcing_snow_depth_for_thf);

   if strcmp(scheme, 'bulk_mo')
      if opts.seb_solver ~= 2
         error('icemodel:configureRun:bulkMoRequiresSebSolver2', ...
            ['turbulent_flux_scheme=''bulk_mo'' currently requires ' ...
            'opts.seb_solver = 2.']);
      end
   end
end

function [vars1, vars2] = defaultOutputVariables(opts)
   % ice1.vars1 = Surface (1-d) data
   % ice2.vars2 = Subsurface (2-d) data
   %
   % Supported output profiles are "standard", "minimal", and "diagnostic".
   %
   % Two important programming notes:
   %
   % 1) The order in which the output variables are set must match the order in
   % which the data are stored in the cell arrays passed to SAVEOUTPUT from the
   % model main functions.
   %
   % 2) If new variables are added, icemodel.postprocess must be reviewed to
   % ensure correct
   % processing is applied, including rounding precision.

   profile = string(opts.output_profile);
   if any(lower(profile) == ["sector", "grid"])
      profile = "minimal";
   end

   switch lower(char(profile))
      case 'minimal'
         if strcmp(opts.smbmodel, 'skinmodel')
            vars1 = {'Tsfc', 'Qm', 'Qe'};
            vars2 = {'Tice'};
         else
            vars1 = {'Tsfc'};
            vars2 = {'Tice', 'f_ice', 'f_liq', 'df_liq', 'df_evp'};
         end

      case 'standard'
         vars1 = ...
            {'Tsfc', 'Qm', 'Qe', 'Qh', 'Qc', 'chi', 'balance', ...
            'dt_sum', 'Tsfc_converged', 'Tice_converged', 'Tice_numiter'};

         if strcmp(opts.smbmodel, 'skinmodel')
            vars2 = {'Tice', 'f_ice', 'f_liq'};
         else
            vars2 = {'Tice', 'f_ice', 'f_liq', 'df_liq', 'df_evp', 'df_lyr', ...
               'Sc', 'r_eff'};
         end

      case 'diagnostic'
         vars1 = ...
            {'Tsfc', 'Qm', 'Qe', 'Qh', 'Qc', 'chi', 'balance', ...
            'dt_sum', 'Tsfc_converged', 'Tice_converged', 'Tice_numiter', ...
            'n_subfail', 'ea_atm', 'De', 'br_coefs_gamma', 'br_coefs_b1_num', ...
            'br_coefs_b2_num', 'roL', 'ro_sfc', ...
            'thf_es_sfc', 'thf_stability_factor', 'thf_z0m', 'thf_z0h', ...
            'thf_z0q', 'thf_u_star', 'thf_L', 'thf_Re', 'thf_numiter'};

         if strcmp(opts.smbmodel, 'skinmodel')
            vars2 = {'Tice', 'f_ice', 'f_liq'};
         else
            vars2 = {'Tice', 'f_ice', 'f_liq', 'df_liq', 'df_evp', 'df_lyr', ...
               'Sc', 'r_eff'};
         end

      otherwise
         error('unrecognized output profile: %s', profile)
   end
end

function opts = configureDebugPaths(opts)
   %CONFIGUREDEBUGPATHS Derive the debug output folder and install env vars.
   %
   % When opts.debug is true and opts.debug_path is empty, the default root
   % is ICEMODEL_OUTPUT_PATH/debug/sitename/smbmodel[/testname], paralleling
   % the standard output path structure. A user-supplied opts.debug_path
   % overrides the entire root.
   %
   % The per-kernel ICEMODEL_DEBUG_*_FILE environment variables are set so
   % the existing dump functions (dumpIceEnbalFailure, dumpMZTransformFailure,
   % etc.) activate without manual env-var configuration.

   if isfield(opts, 'debug_path') && ~isempty(opts.debug_path) ...
         && ~isblanktext(string(opts.debug_path))
      debug_root = opts.debug_path;
   else
      parts = {icemodel.getpath('output'), 'debug', opts.sitename, opts.smbmodel};
      if isfield(opts, 'testname') && ~isempty(opts.testname) ...
            && ~isblanktext(string(opts.testname))
         parts{end+1} = opts.testname;
      end
      debug_root = fullfile(parts{:});
   end

   opts.debug_path = debug_root;

   if exist(debug_root, 'dir') ~= 7
      mkdir(debug_root);
   end

   % Install the per-kernel debug file paths. Each kernel checks its own
   % env var and saves diagnostic state only when the var is non-empty.
   debug_files = { ...
      'ICEMODEL_DEBUG_ICEENBAL_FILE',    'debug_iceenbal.mat'; ...
      'ICEMODEL_DEBUG_MZTRANSFORM_FILE', 'debug_mztransform.mat'; ...
      'ICEMODEL_DEBUG_SEBSOLVE_FILE',    'debug_sebsolve.mat'; ...
      'ICEMODEL_DEBUG_THF_FILE',         'debug_thf.mat'; ...
      'ICEMODEL_DEBUG_SKINSOLVE_FILE',   'debug_skinsolve.mat'; ...
      'ICEMODEL_DEBUG_ICEEBSOLVE_FILE',  'debug_iceebsolve.mat'; ...
      'ICEMODEL_DEBUG_SKINEBSOLVE_FILE', 'debug_skinebsolve.mat'; ...
      'ICEMODEL_DEBUG_MAXSUBSTEP_FILE',  'debug_maxsubstep.mat'; ...
      };
   for k = 1:size(debug_files, 1)
      setenv(debug_files{k, 1}, fullfile(debug_root, debug_files{k, 2}));
   end
end
