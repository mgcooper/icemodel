function opts = configureRun(opts)
   %CONFIGURERUN Fill derived run settings from OPTS.
   %
   %  opts = icemodel.configureRun(opts)
   %
   % This function is idempotent. It fills required derived fields when they
   % are empty, while preserving caller-supplied values for fields such as
   % PATHINPUT, PATHOUTPUT, CASENAME, METFNAME, VARS1, and VARS2.

   if ~isfield(opts, 'pathdata') || isempty(opts.pathdata)
      opts.pathdata = icemodel.setpath('data');
   end

   % OUTPUT_YEARS are the post-spinup years retained in saved/postprocessed
   % output. Keep them derived from SIMYEARS and N_SPINUP_YEARS.
   opts.output_years = icemodel.outputYears(opts);

   if ~isfield(opts, 'pathinput') || isempty(opts.pathinput)
      opts.pathinput = icemodel.setpath('input');
   end
   assert(exist(opts.pathinput, 'dir') == 7, ...
      'ICEMODEL_INPUT_PATH does not exist, set it using icemodel.config');

   if ~isfield(opts, 'patheval') || isempty(opts.patheval)
      opts.patheval = icemodel.setpath('eval');
   end

   if ~isfield(opts, 'pathuserdata') || isempty(opts.pathuserdata)
      opts.pathuserdata = icemodel.setpath('userdata');
   end

   % WRITEOUTPUT appends ['ice1_' opts.casename '.mat'] and saves the file in
   % a subfolder of opts.pathoutput for each year e.g. opts.pathoutput/2016.
   %
   % For grid runs, opts.casename and opts.metfname are set outside this
   % function in the wrapper that loops over grid-cell IDs. Core icemodel only
   % provides the canonical base output folder; any legacy extra subfoldering
   % for gridded workflows belongs in the wrapper, not here.
   if ~isfield(opts, 'pathoutput') || isempty(opts.pathoutput)
      opts.pathoutput = icemodel.setpath('output', ...
         opts.sitename, opts.smbmodel, '', [], opts.testname);
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
end

function [vars1, vars2] = defaultOutputVariables(opts)
   % ice1.vars1 = Surface (1-d) data
   % ice2.vars2 = Subsurface (2-d) data
   %
   % Supported output profiles are "standard" and "minimal".
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
            vars2 = {'Tice', 'f_ice', 'f_liq', 'df_liq', 'df_evp', 'df_lyr', 'Sc'};
         end

      otherwise
         error('unrecognized output profile: %s', profile)
   end
end
