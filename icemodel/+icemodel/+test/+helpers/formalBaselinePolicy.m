function policy = formalBaselinePolicy(baseline_selector)
   %FORMALBASELINEPOLICY Return forcing and default-root policy for a baseline.
   %
   %  policy = icemodel.test.helpers.formalBaselinePolicy("rolling")
   %  policy = icemodel.test.helpers.formalBaselinePolicy("v1.1")
   %  policy = icemodel.test.helpers.formalBaselinePolicy("v1.2")
   %
   % Fields beyond the forcing selection:
   %  promice_filled_policy_sha256  The reconstruction policy digest a frozen
   %     release pins. A blank value means the caller uses the live
   %     icemodel.forcing.reconstruct.policySha256() instead.
   %  use_fixture_root_for_model  True when the release runs the model from
   %     the same provisioned data root it verifies.
   %  require_source_revision  True when the baseline files must record the
   %     source revision that produced them.

   arguments
      baseline_selector (1, :) string = "rolling"
   end

   [baseline_type, baseline_tag] = ...
      icemodel.test.helpers.resolveBaselineSelector(baseline_selector);

   policy = struct();
   policy.baseline_type = baseline_type;
   policy.baseline_tag = baseline_tag;
   % Blank is the default because most baselines have no frozen digest and
   % fall back to the live policy. Each release branch below overrides it.
   policy.promice_filled_policy_sha256 = "";

   % Rolling suites use the verification data tree and one forcing product
   % for every formal station.
   if baseline_type == "rolling"
      policy.config_case = "verification";
      policy.forcing_mode = "fixed";
      policy.forcing = "promice_filled";
      policy.sites = strings(0, 1);
      policy.site_forcings = strings(0, 1);
      policy.required_fixture_capabilities = strings(0, 1);
      policy.use_fixture_root_for_model = false;
      policy.snapshot_from_rolling = false;
      policy.require_source_revision = true;
      policy.promice_filled_policy_sha256 = ...
         icemodel.forcing.reconstruct.policySha256();
      return
   end

   % Preserve the data tree and forcing identity that produced each release's
   % frozen accepted rows.
   release_tag = lower(icemodel.test.helpers.sanitizeTag(baseline_tag));
   switch release_tag
      case "v1_1"
         policy.baseline_tag = "v1.1";
         policy.config_case = "test";
         policy.forcing_mode = "site";
         policy.forcing = "";
         policy.sites = ["kanm"; "kanl"];
         policy.site_forcings = ["kanm"; "kanl"];
         policy.required_fixture_capabilities = "formal-core";
         policy.use_fixture_root_for_model = false;
         policy.snapshot_from_rolling = false;
         policy.require_source_revision = false;

      case "v1_2"
         policy.baseline_tag = "v1.2";
         % use_fixture_root_for_model below makes resolveReleaseDataRoots set
         % data_root to the provisioned test/data tree, which overrides this
         % config case for every suite entry point.
         policy.config_case = "verification";
         policy.forcing_mode = "fixed";
         policy.forcing = "promice_filled";
         policy.sites = strings(0, 1);
         policy.site_forcings = strings(0, 1);
         policy.required_fixture_capabilities = "formal-core";
         policy.use_fixture_root_for_model = true;
         policy.snapshot_from_rolling = true;
         policy.require_source_revision = true;
         policy.promice_filled_policy_sha256 = ...
            "bd336da0880474f1987facc2311c4f45a6c281877ae8b944a3fbdc7cfb68d513";

      otherwise
         error('icemodel:test:unregisteredReleaseForcing', ...
            ['Release baseline %s has no registered formal forcing identity. ', ...
            'Register it in formalBaselinePolicy before comparison or build.'], ...
            baseline_tag)
   end
end
