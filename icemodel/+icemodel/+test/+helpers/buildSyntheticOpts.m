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
      kwargs.seb_solver (1, 1) double = NaN
      kwargs.turbulent_flux_scheme (1, :) char = ''
      kwargs.z0_bulk (1, 1) double = NaN
      kwargs.z0_ice (1, 1) double = NaN
      kwargs.z0_snow_low_density (1, 1) double = NaN
      kwargs.z0_snow_high_density (1, 1) double = NaN
      kwargs.n_spinup_years (1, 1) double {mustBeNonnegative, mustBeInteger} = 0
      kwargs.pathoutput (1, :) char = ''
      kwargs.use_restart (1, 1) logical = false
      kwargs.restartfile (1, :) char = ''
      kwargs.saverestart (1, 1) logical = false
      kwargs.metfname = {}
   end

   % Test fixtures should be explicit about the THF scheme they exercise.
   % Only tests that intentionally check the runtime default should depend on
   % the ambient setopts default.
   turbulent_flux_scheme = kwargs.turbulent_flux_scheme;
   if isempty(turbulent_flux_scheme)
      turbulent_flux_scheme = 'bulk_richardson';
   end

   % Start from the normal model defaults for this synthetic site/model, but
   % pin the THF scheme up front so setopts validation never depends on an
   % ambient user-edited default.
   setopts_args = {'turbulent_flux_scheme', turbulent_flux_scheme};
   if strcmpi(turbulent_flux_scheme, 'bulk_mo')
      if isfinite(kwargs.solver)
         setopts_args(end+1:end+2) = {'solver', kwargs.solver};
      else
         setopts_args(end+1:end+2) = {'solver', 1};
      end
      if isfinite(kwargs.seb_solver)
         setopts_args(end+1:end+2) = {'seb_solver', kwargs.seb_solver};
      else
         setopts_args(end+1:end+2) = {'seb_solver', 2};
      end
   end
   opts = icemodel.setopts(smbmodel, workspace.sitename, simyears, ...
      workspace.forcings, kwargs.userdata, kwargs.uservars, ...
      kwargs.testname, kwargs.saveflag, kwargs.backupflag, setopts_args{:});

   % Prefer the fixture timestep unless the caller overrode it explicitly.
   dt_value = kwargs.dt;
   if ~isfinite(dt_value)
      dt_value = workspace.dt_seconds;
   end

   % Collect the override fields once so RESETOPTS stays data-driven.
   names = {'n_spinup_years', 'pathinput', 'pathuserdata', 'patheval', ...
      'output_profile', 'use_restart', 'restartfile', 'saverestart', 'dt', ...
      'turbulent_flux_scheme'};
   values = {kwargs.n_spinup_years, workspace.inputdir, ...
      workspace.userdatadir, workspace.evaldir, kwargs.output_profile, ...
      kwargs.use_restart, kwargs.restartfile, kwargs.saverestart, dt_value, ...
      turbulent_flux_scheme};
   if isfinite(kwargs.solver)
      names{end+1} = 'solver';
      values{end+1} = kwargs.solver;
   end
   if isfinite(kwargs.seb_solver)
      names{end+1} = 'seb_solver';
      values{end+1} = kwargs.seb_solver;
   end
   if isfinite(kwargs.z0_bulk)
      names{end+1} = 'z0_bulk';
      values{end+1} = kwargs.z0_bulk;
   end
   if isfinite(kwargs.z0_ice)
      names{end+1} = 'z0_ice';
      values{end+1} = kwargs.z0_ice;
   end
   if isfinite(kwargs.z0_snow_low_density)
      names{end+1} = 'z0_snow_low_density';
      values{end+1} = kwargs.z0_snow_low_density;
   end
   if isfinite(kwargs.z0_snow_high_density)
      names{end+1} = 'z0_snow_high_density';
      values{end+1} = kwargs.z0_snow_high_density;
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
