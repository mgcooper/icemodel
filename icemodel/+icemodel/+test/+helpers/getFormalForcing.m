function forcing = getFormalForcing(kwargs)
   %GETFORMALFORCING Return the forcing identity for one formal baseline.
   %
   %  forcing = icemodel.test.helpers.getFormalForcing()
   %  forcing = icemodel.test.helpers.getFormalForcing( ...
   %     sitename="kanm", baseline="v1.1")
   %
   % Rolling formal suites use the official gap-filled PROMICE forcing.
   % Frozen release comparisons retain the forcing identity with which their
   % accepted rows were produced. Release tags must be registered in
   % formalBaselinePolicy before their case matrices can run, so a new tag
   % cannot silently inherit an unrelated forcing product.

   arguments
      kwargs.sitename (1, 1) string = ""
      kwargs.baseline (1, :) string = "rolling"
   end

   policy = icemodel.test.helpers.formalBaselinePolicy(kwargs.baseline);

   % The mutable/default suite follows the accepted forcing going forward.
   if policy.forcing_mode == "fixed"
      forcing = policy.forcing;
      return
   end

   % Frozen v1.1 rows predate promice_filled and used their station-specific
   % forcing aliases. Keep the supported aliases explicit rather than treating
   % an arbitrary site name as a runnable forcing.
   idx = find(strcmpi(kwargs.sitename, policy.sites), 1);
   if isempty(idx)
      error('icemodel:test:unsupportedReleaseForcingSite', ...
         'Release baseline %s has no registered forcing for site %s.', ...
         policy.baseline_tag, kwargs.sitename)
   end
   forcing = policy.site_forcings(idx);
end
