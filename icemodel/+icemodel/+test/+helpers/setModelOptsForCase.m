function opts = setModelOptsForCase(c, kwargs)
   %SETMODELOPTSFORCASE Build resolved model OPTS for one formal test case.
   %
   %  opts = icemodel.test.helpers.setModelOptsForCase(c)
   %  opts = icemodel.test.helpers.setModelOptsForCase(c, include_spinup=false)
   %
   % Uses the requested model defaults for the case, then applies the explicit
   % case-specific reset(s) needed by the formal matrix.
   arguments
      c
      kwargs.include_spinup (1, 1) logical = true
   end

   % Accept either a one-row case table or an equivalent scalar struct.
   if istable(c)
      assert(height(c) == 1, ...
         'setModelOptsForCase expects a single table row')
      c = table2struct(c);
   end

   % Expand the case into explicit simulation years and spinup policy.
   [simyears, n_spinup_years] = resolveSimulationYears(c, kwargs.include_spinup);

   % Start from the public runtime contract before applying case-specific
   % overrides from the formal matrix.
   opts = icemodel.setopts(c.smbmodel, c.sitename, simyears, ...
      c.forcings, c.userdata, c.uservars, string.empty(), false, false, ...
      'n_spinup_years', n_spinup_years);

   % Solver selection is the main explicit override used by the suites.
   if isfield(c, 'solver') && ~isempty(c.solver)
      opts = icemodel.resetopts(opts, 'solver', c.solver);
   end

end

function [simyears, n_spinup_years] = resolveSimulationYears(c, include_spinup)
   %RESOLVESIMULATIONYEARS Resolve retained and spinup years for one case.

   % Formal cases may carry either one SIMYEAR or an explicit SIMYEARS
   % vector for spinup-aware comparisons.
   if isfield(c, 'simyears') && ~isempty(c.simyears)
      simyears = c.simyears;
   elseif isfield(c, 'simyear') && ~isempty(c.simyear)
      simyears = [c.simyear - 1, c.simyear];
   else
      error('formal case is missing simyear/simyears')
   end

   % Default formal cases to one leading spinup year when only SIMYEAR is
   % given and the case matrix has not overridden that policy explicitly.
   if isfield(c, 'n_spinup_years') && ~isempty(c.n_spinup_years)
      n_spinup_years = c.n_spinup_years;
   elseif isfield(c, 'simyear') && ~isempty(c.simyear)
      n_spinup_years = 1;
   else
      n_spinup_years = 0;
   end

   % Perf-style callers can drop the spinup years from the runtime contract
   % while still reusing the same formal case definition.
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
