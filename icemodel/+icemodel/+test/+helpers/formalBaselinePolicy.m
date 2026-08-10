function policy = formalBaselinePolicy(baseline_selector)
   %FORMALBASELINEPOLICY Return forcing and default-root policy for a baseline.
   %
   %  policy = icemodel.test.helpers.formalBaselinePolicy("rolling")
   %  policy = icemodel.test.helpers.formalBaselinePolicy("v1.1")

   arguments
      baseline_selector (1, :) string = "rolling"
   end

   [baseline_type, baseline_tag] = ...
      icemodel.test.helpers.resolveBaselineSelector(baseline_selector);

   policy = struct();
   policy.baseline_type = baseline_type;
   policy.baseline_tag = baseline_tag;

   % Rolling suites use the verification data tree and one forcing product
   % for every formal station.
   if baseline_type == "rolling"
      policy.config_case = "verification";
      policy.forcing_mode = "fixed";
      policy.forcing = "promice_filled";
      policy.sites = strings(0, 1);
      policy.site_forcings = strings(0, 1);
      policy.required_fixture_capabilities = strings(0, 1);
      policy.snapshot_from_rolling = false;
      return
   end

   % Release registrations preserve both the data tree and forcing identity
   % used to produce their frozen accepted rows.
   release_tag = lower(icemodel.test.helpers.sanitizeTag(baseline_tag));
   switch release_tag
      case "v1_1"
         policy.config_case = "test";
         policy.forcing_mode = "site";
         policy.forcing = "";
         policy.sites = ["kanm"; "kanl"];
         policy.site_forcings = ["kanm"; "kanl"];
         policy.required_fixture_capabilities = "formal-core";
         policy.snapshot_from_rolling = false;

      otherwise
         error('icemodel:test:unregisteredReleaseForcing', ...
            ['Release baseline %s has no registered formal forcing identity. ', ...
            'Register it in formalBaselinePolicy before comparison or build.'], ...
            baseline_tag)
   end
end
