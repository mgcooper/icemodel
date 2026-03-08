function opts = configureRun(opts)
   %CONFIGURERUN Fill derived run settings from OPTS.
   %
   %  opts = icemodel.configureRun(opts)
   %
   % This function is idempotent. It fills required derived fields when they
   % are empty, while preserving caller-supplied values for fields such as
   % PATHINPUT, PATHOUTPUT, CASENAME, METFNAME, VARS1, and VARS2.

   if ~isfield(opts, 'pathinput') || isempty(opts.pathinput)
      opts.pathinput = getenv('ICEMODEL_INPUT_PATH');
   end
   assert(exist(opts.pathinput, 'dir') == 7, ...
      'ICEMODEL_INPUT_PATH does not exist, set it using icemodel.config');

   % WRITEOUTPUT appends ['ice1_' opts.casename '.mat'] and saves the file in
   % a subfolder of opts.pathoutput for each year e.g. opts.pathoutput/2016.
   %
   % For gridded sector runs, opts.casename is set to the grid point ID
   % outside this function in a loop (icemodel.run.grid in the runoff project),
   % to get ice1_1.mat, ice1_2.mat, and so forth. Same for metfname. Writing to
   % netcdf should eliminate this.
   %
   % Also note that for gridded sector runs, "userdata" is appended to the
   % output folder. This is a legacy method: sector output filenames use the
   % grid-cell number rather than the forcing/userdata/uservars naming used for
   % point runs, so the swap case is distinguished at the folder level:
   %    ICEMODEL_OUTPUT_PATH/sector/icemodel/mar
   %    ICEMODEL_OUTPUT_PATH/sector/icemodel/modis
   % In practice these sector runs used MAR forcing in both cases, with the
   % second folder indicating the modis-albedo swap.
   if ~isfield(opts, 'pathoutput') || isempty(opts.pathoutput)
      opts.pathoutput = fullfile( ...
         getenv('ICEMODEL_OUTPUT_PATH'), opts.sitename, opts.smbmodel);
      if strcmp(opts.sitename, 'sector')
         opts.pathoutput = fullfile(opts.pathoutput, opts.userdata);
      end
      opts.pathoutput = fullfile(opts.pathoutput, opts.testname);
   end

   % Create the casename. WRITEOUTPUT appends this to the base filenames.
   % Note that for gridded runs, casename is modified to represent the grid
   % number, see icemodel.run.grid in the runoff project.
   if ~isfield(opts, 'casename') || isempty(opts.casename)
      opts.casename = icemodel.setcase(opts.forcings, opts.userdata, opts.uservars);
   end

   if ~isfield(opts, 'output_profile') || isempty(opts.output_profile)
      if strcmp(opts.sitename, 'sector')
         opts.output_profile = 'minimal';
      else
         opts.output_profile = 'standard';
      end
   end

   if ~isfield(opts, 'metfname') || isempty(opts.metfname)
      if strcmp(opts.sitename, 'sector')
         % Could add logic here to deal with sector file names. For now, the
         % metfname must be set outside this function in a loop.
         % for n = 1:numel(runpoints)
         %    opts.metfname = 'met_sector.mat';
         % end
         opts.metfname = {};
      else
         opts.metfname = fullfile(opts.pathinput, 'met', ...
            icemodel.createMetFileNames(opts));
      end
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
   % 2) If new variables are added, POSTPROC must be reviewed to ensure correct
   % processing is applied, including rounding precision.

   profile = string(opts.output_profile);
   if lower(profile) == "sector"
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
