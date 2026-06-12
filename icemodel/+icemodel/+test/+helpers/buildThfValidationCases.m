function cases = buildThfValidationCases(kwargs)
   %BUILDTHFVALIDATIONCASES Build focused real-case THF validation cases.
   %
   %  cases = icemodel.test.helpers.buildThfValidationCases()
   %  cases = icemodel.test.helpers.buildThfValidationCases(...)
   %
   % Name-value options:
   %   sites        - string array of site names, default ["kanm" "kanl"]
   %   schemes      - string array of THF schemes, default
   %                  ["bulk_richardson" "monin_obukhov"]
   %   spinup_sites - subset of SITES that should also include the
   %                  [2015 2016] / n_spinup_years = 1 cases

   arguments
      kwargs.sites string = ["kanm" "kanl"]
      kwargs.schemes string = ["bulk_richardson" "monin_obukhov"]
      kwargs.spinup_sites string = string.empty(1, 0)
   end

   sites = lower(string(kwargs.sites));
   schemes = lower(string(kwargs.schemes));
   spinup_sites = lower(string(kwargs.spinup_sites));

   % Upper bound: each site contributes one case per scheme, plus an
   % optional spinup case per scheme when the site is in spinup_sites.
   max_cases = numel(sites) * numel(schemes) * 2;
   cases = repmat(initCase(), max_cases, 1);
   n = 0;
   for i_site = 1:numel(sites)
      site = sites(i_site);

      for i_scheme = 1:numel(schemes)
         scheme = schemes(i_scheme);
         n = n + 1;
         cases(n, 1) = makeCase(site, 2016, 0, scheme);
      end

      if any(spinup_sites == site)
         for i_scheme = 1:numel(schemes)
            scheme = schemes(i_scheme);
            n = n + 1;
            cases(n, 1) = makeCase(site, [2015 2016], 1, scheme);
         end
      end
   end
   cases = cases(1:n, 1);
end

function c = initCase()
   %INITCASE Initialize one focused THF validation case struct.

   c = struct( ...
      'label', "", ...
      'site', "", ...
      'forcings', "", ...
      'simyears', [], ...
      'n_spinup_years', 0, ...
      'scheme', "");
end

function c = makeCase(site, simyears, n_spinup_years, scheme)
   %MAKECASE Build one focused THF validation case struct.

   c = initCase();
   c.site = site;
   c.forcings = site;
   c.simyears = simyears;
   c.n_spinup_years = n_spinup_years;
   c.scheme = scheme;
   c.label = composeLabel(site, simyears, n_spinup_years, scheme);
end

function label = composeLabel(site, simyears, n_spinup_years, scheme)
   %COMPOSELABEL Build the canonical human-readable case label.

   years_label = join(string(simyears), "-");
   site_label = upper(char(site));

   if n_spinup_years > 0
      label = string(sprintf('%s %s spinup %s', ...
         site_label, years_label, scheme));
   else
      label = string(sprintf('%s %s %s', ...
         site_label, years_label, scheme));
   end
end
