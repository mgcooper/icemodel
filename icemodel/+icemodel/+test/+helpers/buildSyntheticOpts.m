function opts = buildSyntheticOpts(workspace, smbmodel, simyears, kwargs)
   %BUILDSYNTHETICOPTS Build resolved OPTS for synthetic unit-test runs.
   %
   %  opts = icemodel.test.helpers.buildSyntheticOpts(workspace, "skinmodel", 2016)

   arguments
      workspace struct
      smbmodel (1, :) char {mustBeMember(smbmodel, {'icemodel', 'skinmodel'})}
      simyears
      kwargs.userdata (1, :) char = ''
      kwargs.uservars (1, :) char = ''
      kwargs.testname (1, :) char = ''
      kwargs.saveflag (1, 1) logical = false
      kwargs.backupflag (1, 1) logical = false
      kwargs.output_profile (1, :) char = 'standard'
      kwargs.dt (1, 1) double = NaN
      kwargs.solver (1, 1) double = NaN
      kwargs.n_spinup_years (1, 1) double {mustBeNonnegative, mustBeInteger} = 0
      kwargs.pathoutput (1, :) char = ''
      kwargs.use_restart (1, 1) logical = false
      kwargs.restartfile (1, :) char = ''
      kwargs.saverestart (1, 1) logical = false
      kwargs.metfname = {}
   end

   % Start from the normal model defaults for this synthetic site/model.
   opts = icemodel.setopts(smbmodel, workspace.sitename, simyears, ...
      workspace.forcings, kwargs.userdata, kwargs.uservars, ...
      kwargs.testname, kwargs.saveflag, kwargs.backupflag);

   % Prefer the fixture timestep unless the caller overrode it explicitly.
   dt_value = kwargs.dt;
   if ~isfinite(dt_value)
      dt_value = workspace.dt_seconds;
   end

   % Collect the override fields once so RESETOPTS stays data-driven.
   names = {'n_spinup_years', 'pathinput', 'pathuserdata', 'patheval', ...
      'output_profile', 'use_restart', 'restartfile', 'saverestart', 'dt'};
   values = {kwargs.n_spinup_years, workspace.inputdir, ...
      workspace.userdatadir, workspace.evaldir, kwargs.output_profile, ...
      kwargs.use_restart, kwargs.restartfile, kwargs.saverestart, dt_value};
   if isfinite(kwargs.solver)
      names{end+1} = 'solver';
      values{end+1} = kwargs.solver;
   end
   if ~isempty(kwargs.pathoutput)
      names{end+1} = 'pathoutput';
      values{end+1} = kwargs.pathoutput;
   end
   if ~isempty(kwargs.metfname)
      names{end+1} = 'metfname';
      values{end+1} = kwargs.metfname;
   end

   % Finalize all derived run fields before returning the OPTS struct.
   resetargs = reshape([names; values], 1, []);
   opts = icemodel.resetopts(opts, resetargs{:});
   opts = icemodel.configureRun(opts);
end
