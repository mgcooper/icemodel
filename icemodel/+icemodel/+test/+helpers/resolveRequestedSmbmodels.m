function models = resolveRequestedSmbmodels(smbmodel)
   %RESOLVEREQUESTEDSMBMODELS Expand one requested formal smbmodel selector.
   %
   %  models = icemodel.test.helpers.resolveRequestedSmbmodels("all")
   %  models = icemodel.test.helpers.resolveRequestedSmbmodels("icemodel")
   %
   % Use this helper at a formal-suite entry point. The workflow then runs one
   % named model, and loops over every model when the caller asks for the
   % aggregate selector "all".

   arguments
      smbmodel (1, :) string ...
         {icemodel.validators.mustBeTestSmbmodelSelector(smbmodel)}
   end

   % Expand the aggregate selector once, at the entry point.
   if smbmodel == "all"
      models = icemodel.namelists.smbmodel("test");
   else
      models = reshape(smbmodel, [], 1);
   end

   % Return a column vector so callers can loop without shape checks.
   models = string(models(:));
end
