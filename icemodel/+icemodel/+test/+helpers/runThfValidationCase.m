function [ice1, met, opts, metrics, runtime_seconds] = runThfValidationCase(c, kwargs)
   %RUNTHFVALIDATIONCASE Run one real-case THF validation scenario.
   %
   %  [ice1, met, opts, metrics, runtime_seconds] = ...
   %     icemodel.test.helpers.runThfValidationCase(c)
   %  [...] = icemodel.test.helpers.runThfValidationCase(c, overrides=S)
   %
   % Input C is a scalar struct or one-row table describing the validation
   % case. The required fields are:
   %   site or sitename      - site identifier such as "kanm" or "kanl"
   %   simyears or simyear   - retained years or explicit spinup+retained years
   %
   % Optional fields are:
   %   forcings             - forcing source; defaults to site/sitename
   %   smbmodel             - defaults to "icemodel"
   %   n_spinup_years       - defaults to 0
   %   scheme               - defaults to "bulk_richardson"
   %
   % OVERRIDES is a struct of explicit OPTS overrides, such as:
   %   struct('z0_ice', 0.001, 'testname', 'kanl_bulk_mo_rough')
   %
   % This helper applies the accepted THF validation contract by default:
   % `solver = 1` and `seb_solver = 2`. Callers can still override those
   % fields explicitly via OVERRIDES when a validation pass needs to probe a
   % different runtime configuration.

   arguments
      c
      kwargs.overrides struct = struct()
   end

   if istable(c)
      assert(height(c) == 1, ...
         'runThfValidationCase expects a scalar struct or one-row table')
      c = table2struct(c);
   end

   sitename = getCaseField(c, {'site', 'sitename'});
   if isfield(c, 'forcings') && ~isempty(c.forcings)
      forcings = c.forcings;
   else
      forcings = sitename;
   end

   if isfield(c, 'simyears') && ~isempty(c.simyears)
      simyears = c.simyears;
   else
      simyears = getCaseField(c, {'simyear'});
   end

   if isfield(c, 'smbmodel') && ~isempty(c.smbmodel)
      smbmodel = c.smbmodel;
   else
      smbmodel = 'icemodel';
   end

   if isfield(c, 'n_spinup_years') && ~isempty(c.n_spinup_years)
      n_spinup_years = c.n_spinup_years;
   else
      n_spinup_years = 0;
   end

   if isfield(c, 'scheme') && ~isempty(c.scheme)
      scheme = c.scheme;
   else
      scheme = 'bulk_richardson';
   end

   opts = icemodel.setopts(smbmodel, sitename, simyears, forcings, ...
      [], [], [], false, false, ...
      'n_spinup_years', n_spinup_years, ...
      'solver', 1, ...
      'seb_solver', 2, ...
      'turbulent_flux_scheme', scheme);

   opts = applyOverrides(opts, kwargs.overrides);

   tic
   [ice1, ice2, opts] = icemodel(opts);
   runtime_seconds = toc;

   [ice1, ~, met] = icemodel.postprocess(ice1, ice2, opts, opts.output_years);
   metrics = icemodel.test.helpers.summarizeIce1Metrics(ice1, met);
end

function value = getCaseField(c, names)
   %GETCASEFIELD Return the first matching field from the case struct.

   for i = 1:numel(names)
      name = names{i};
      if isfield(c, name) && ~isempty(c.(name))
         value = c.(name);
         return
      end
   end

   error('runThfValidationCase:missingField', ...
      'Missing required case field: %s', strjoin(names, ' or '));
end

function opts = applyOverrides(opts, overrides)
   %APPLYOVERRIDES Apply explicit resetopts overrides from one struct.

   names = fieldnames(overrides);
   if isempty(names)
      return
   end

   args = cell(1, 2 * numel(names));
   for i = 1:numel(names)
      args{2*i - 1} = names{i};
      args{2*i} = overrides.(names{i});
   end
   opts = icemodel.resetopts(opts, args{:});
end
