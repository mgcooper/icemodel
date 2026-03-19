function models = resolveRequestedSmbmodels(smbmodel)
   %RESOLVEREQUESTEDSMBMODELS Expand one requested formal smbmodel selector.
   %
   %  models = icemodel.test.helpers.resolveRequestedSmbmodels("all")
   %  models = icemodel.test.helpers.resolveRequestedSmbmodels("icemodel")
   %
   % Use this helper at formal-suite entrypoints so the canonical workflow is
   % always "run one concrete model, then loop when the caller requested the
   % virtual aggregate selector".

   arguments
      smbmodel (1, :) string ...
         {icemodel.validators.mustBeTestSmbmodelSelector(smbmodel)}
   end

   % Expand the virtual aggregate selector once at the entrypoint boundary.
   if smbmodel == "all"
      models = icemodel.namelists.smbmodel("test");
   else
      models = reshape(smbmodel, [], 1);
   end

   % Return a column vector so callers can loop without shape checks.
   models = string(models(:));
end
