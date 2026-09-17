function [data_root, fixture_root] = resolveReleaseDataRoots( ...
      policy, data_root, fixture_root)
   %RESOLVERELEASEDATAROOTS Resolve model and fixture roots for a baseline.
   %
   %  [data_root, fixture_root] = ...
   %     icemodel.test.helpers.resolveReleaseDataRoots( ...
   %     policy, data_root, fixture_root)
   %
   % See also: icemodel.test.helpers.formalBaselinePolicy,
   %  icemodel.verification.setup.fixtureDataRoot

   arguments
      policy (1, 1) struct
      data_root (1, 1) string
      fixture_root (1, 1) string
   end

   % A baseline that provisions nothing has no fixture root to resolve, so
   % both caller values stand.
   if isempty(policy.required_fixture_capabilities)
      return
   end

   % Verify inside the tree the caller named. Fall back to the release's own
   % provisioned root only when the caller named no tree at all.
   if isblanktext(fixture_root)
      if ~isblanktext(data_root)
         fixture_root = data_root;
      else
         fixture_root = icemodel.verification.setup.fixtureDataRoot( ...
            policy.baseline_tag);
      end
   end

   % A release that runs the model from its provisioned data must verify and
   % execute the same tree. Otherwise the comparison would gate accepted rows
   % on data the hash check never saw.
   if policy.use_fixture_root_for_model
      if isblanktext(data_root)
         data_root = fixture_root;
      elseif icemodel.helpers.canonicalPath(data_root) ...
            ~= icemodel.helpers.canonicalPath(fixture_root)
         error('icemodel:test:releaseDataRootMismatch', ...
            ['This release must verify and execute the same ' ...
            'provisioned data root.'])
      end
   end
end
