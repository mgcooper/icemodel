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

   if istable(c)
      assert(height(c) == 1, ...
         'setModelOptsForCase expects a single table row')
      c = table2struct(c);
   end

   [simyears, n_spinup_years] = resolveSimulationYears(c, kwargs.include_spinup);

   opts = icemodel.setopts(c.smbmodel, c.sitename, simyears, ...
      c.forcings, c.userdata, c.uservars, string.empty(), false, false, ...
      'n_spinup_years', n_spinup_years);

   if isfield(c, 'solver') && ~isempty(c.solver)
      opts = icemodel.resetopts(opts, 'solver', c.solver);
   end
end

function [simyears, n_spinup_years] = resolveSimulationYears(c, include_spinup)

   if isfield(c, 'simyears') && ~isempty(c.simyears)
      simyears = c.simyears;
   elseif isfield(c, 'simyear') && ~isempty(c.simyear)
      simyears = c.simyear;
   else
      error('formal case is missing simyear/simyears')
   end

   if isfield(c, 'n_spinup_years') && ~isempty(c.n_spinup_years)
      n_spinup_years = c.n_spinup_years;
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
