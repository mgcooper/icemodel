function opts = buildFormalCaseOpts(c)
%BUILDFORMALCASEOPTS Build resolved OPTS for one formal icemodel test case.
%
%  opts = test.helpers.buildFormalCaseOpts(c)
%
% Uses the requested model defaults for the case, then applies the explicit
% case-specific reset(s) needed by the formal matrix.
   arguments
      c
   end

   if istable(c)
      assert(height(c) == 1, ...
         'buildFormalCaseOpts expects a single table row')
      c = table2struct(c);
   end

   opts = icemodel.setopts(c.smbmodel, c.sitename, c.simyear, ...
      c.forcings, c.userdata, c.uservars, string.empty(), false, false);

   if isfield(c, 'solver_mode') && ~isempty(c.solver_mode)
      opts = icemodel.resetopts(opts, 'bc_type', c.solver_mode);
   end
end
