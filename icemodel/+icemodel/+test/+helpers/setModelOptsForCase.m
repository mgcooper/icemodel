function opts = setModelOptsForCase(c, kwargs)
   %SETMODELOPTSFORCASE Build resolved model OPTS for one case.
   %
   %  opts = icemodel.test.helpers.setModelOptsForCase(c)
   %  opts = icemodel.test.helpers.setModelOptsForCase(c, include_spinup=false)
   %  opts = icemodel.test.helpers.setModelOptsForCase(case_manifest)
   %  opts = icemodel.test.helpers.setModelOptsForCase(case_manifest, ...
   %     startdate=..., enddate=..., dt=900)
   %
   % Accepts either shape and dispatches:
   %
   %   Formal case (regression / perf / baseline matrices)
   %     Row table or scalar struct with fields smbmodel, sitename, forcings,
   %     userdata, uservars, simyear or simyears, optional n_spinup_years
   %     and solver. Resolves simyears + spinup policy, then applies
   %     case-specific overrides on top of the public icemodel.setopts
   %     contract.
   %
   %   Verification manifest (verification suite)
   %     Scalar struct with case_id and comparison_window. Uses case_id as
   %     both sitename and forcings, derives simyears from the comparison
   %     window, and carries the window through opts.startdate / opts.enddate
   %     so the standard configureRun + createMetFileNames + loadmet chain
   %     resolves the staged multi-year met file with no special-case bypass.
   %
   % Inputs
   %   c   Formal case row/struct OR verification case manifest.
   %
   % Name-value
   %   include_spinup : logical, default true. Formal-case only. When false,
   %       drop the n_spinup_years leading years from simyears.
   %   smbmodel : string, default "icemodel". Verification-manifest only.
   %   testname : string, default "verification". Verification-manifest only.
   %   dt : double seconds. Verification-manifest only. Default derived from
   %       manifest.native_timestep ("hourly" -> 3600).
   %   startdate / enddate : datetime, default NaT. Verification-manifest only.
   %       Explicit window narrowing; defaults to manifest.comparison_window.
   %
   % See also: icemodel.setopts, icemodel.createMetFileNames,
   %  icemodel.verification.runIcemodelSnowCandidate

   arguments
      c
      kwargs.include_spinup (1, 1) logical = true
      kwargs.smbmodel (1, 1) string = "icemodel"
      kwargs.testname (1, 1) string = "verification"
      kwargs.dt double = []
      kwargs.startdate = NaT('TimeZone', 'UTC')
      kwargs.enddate = NaT('TimeZone', 'UTC')
   end

   % Accept either a one-row case table or an equivalent scalar struct.
   if istable(c)
      assert(height(c) == 1, ...
         'setModelOptsForCase expects a single table row')
      c = table2struct(c);
   end

   % Dispatch on input shape. Verification manifests carry case_id and a
   % comparison_window; formal cases carry simyear / simyears plus the
   % full positional contract for icemodel.setopts.
   if isfield(c, 'case_id') && isfield(c, 'comparison_window')
      opts = optsFromVerificationManifest(c, kwargs);
   else
      opts = optsFromFormalCase(c, kwargs);
   end
end

function opts = optsFromFormalCase(c, kwargs)
   %OPTSFROMFORMALCASE Build opts from a formal regression/perf case row.

   [simyears, n_spinup_years] = resolveSimulationYears(c, kwargs.include_spinup);

   opts = icemodel.setopts(c.smbmodel, c.sitename, simyears, ...
      c.forcings, c.userdata, c.uservars, string.empty(), false, false, ...
      'n_spinup_years', n_spinup_years);

   if isfield(c, 'solver') && ~isempty(c.solver)
      opts = icemodel.resetopts(opts, 'solver', c.solver);
   end
end

function opts = optsFromVerificationManifest(case_manifest, kwargs)
   %OPTSFROMVERIFICATIONMANIFEST Build opts from a verification case manifest.

   [window_start, window_end] = resolveComparisonWindow(case_manifest, ...
      kwargs.startdate, kwargs.enddate);

   dt_seconds = resolveTimestep(case_manifest, kwargs.dt);

   smbmodel = char(kwargs.smbmodel);
   sitename = char(string(case_manifest.case_id));
   simyears = unique(year([window_start; window_end]))';
   forcings = sitename;
   userdata = [];
   uservars = [];
   testname = char(kwargs.testname);
   saveflag = false;
   backupflag = false;

   opts = icemodel.setopts(smbmodel, sitename, simyears, forcings, ...
      userdata, uservars, testname, saveflag, backupflag, ...
      'dt', dt_seconds, ...
      'startdate', window_start, ...
      'enddate', window_end);
end

function [simyears, n_spinup_years] = resolveSimulationYears(c, include_spinup)
   %RESOLVESIMULATIONYEARS Resolve retained and spinup years for one formal case.

   if isfield(c, 'simyears') && ~isempty(c.simyears)
      simyears = c.simyears;
   elseif isfield(c, 'simyear') && ~isempty(c.simyear)
      simyears = [c.simyear - 1, c.simyear];
   else
      error('formal case is missing simyear/simyears')
   end

   if isfield(c, 'n_spinup_years') && ~isempty(c.n_spinup_years)
      n_spinup_years = c.n_spinup_years;
   elseif isfield(c, 'simyear') && ~isempty(c.simyear)
      n_spinup_years = 1;
   else
      n_spinup_years = 0;
   end

   if ~include_spinup
      if n_spinup_years >= numel(simyears)
         error('n_spinup_years must be smaller than numel(simyears)')
      end
      if n_spinup_years > 0
         simyears = simyears(n_spinup_years+1:end);
      end
      n_spinup_years = 0;
   end
end

function [window_start, window_end] = resolveComparisonWindow(manifest, ...
      startdate, enddate)
   %RESOLVECOMPARISONWINDOW Resolve the canonical window for a manifest case.

   if ~isnat(startdate) && ~isnat(enddate)
      window_start = icemodel.verification.setup.ensureUtc(startdate);
      window_end = icemodel.verification.setup.ensureUtc(enddate);
      return
   end

   window_start = icemodel.verification.setup.ensureUtc( ...
      manifest.comparison_window.start);
   window_end = icemodel.verification.setup.ensureUtc( ...
      manifest.comparison_window.end);
end

function dt_seconds = resolveTimestep(case_manifest, dt_override)
   %RESOLVETIMESTEP Map the manifest's native_timestep label to dt seconds.

   if ~isempty(dt_override)
      dt_seconds = dt_override;
      return
   end
   if isfield(case_manifest, 'native_timestep')
      label = lower(string(case_manifest.native_timestep));
   else
      label = "";
   end
   switch label
      case "hourly"
         dt_seconds = 3600;
      case {"15-min", "15 min", "15m", "15-minute", "15 minute"}
         dt_seconds = 900;
      otherwise
         dt_seconds = 3600;
   end
end
