function validateEvalTarget(eval_target)
   %VALIDATEEVALTARGET Validate a case-manifest eval_target value.
   %
   %  icemodel.verification.setup.validateEvalTarget(eval_target)
   %
   % Inputs
   %  eval_target   The eval_target value from a case-manifest entry. A string
   %                array naming which model capabilities the case exercises. An
   %                empty value (string(0,1) or "") is permitted (cases that
   %                exercise no curated capability, e.g. analytical benchmarks);
   %                every non-empty element must be a member of the canonical
   %                vocabulary published by icemodel.verification.namelists.evaltarget.
   %
   % Role
   %  Setup-side schema gate shared by makeCaseManifestEntry and
   %  makeFirnCaseManifestEntry so a stamped eval_target cannot drift from the
   %  canonical namelist.
   %
   % See also: icemodel.verification.namelists.evaltarget,
   %  icemodel.verification.setup.validateSurfaceZone

   targets = string(eval_target);
   targets = targets(strlength(targets) > 0 & ~ismissing(targets));
   if isempty(targets)
      return
   end

   allowed = icemodel.verification.namelists.evaltarget();
   bad = targets(~ismember(targets, allowed));
   if ~isempty(bad)
      error('icemodel:verification:setup:invalidEvalTarget', ...
         'eval_target "%s" is not in the canonical vocabulary (%s)', ...
         strjoin(unique(bad), ', '), strjoin(allowed, ', '))
   end
end
