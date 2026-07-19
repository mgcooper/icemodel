function validatePermafrostZone(permafrost_zone)
   %VALIDATEPERMAFROSTZONE Validate a case-manifest permafrost_zone value.
   %
   %  icemodel.verification.setup.validatePermafrostZone(permafrost_zone)
   %
   % Inputs
   %  permafrost_zone   The permafrost_zone value from a case-manifest entry. An
   %                    empty value ("") is permitted (cases where the permafrost
   %                    extent is not meaningful); any non-empty value must be a
   %                    member of the canonical vocabulary published by
   %                    icemodel.verification.namelists.permafrostzone.
   %
   % Role
   %  Setup-side schema gate shared by makeCaseManifestEntry and
   %  makeFirnCaseManifestEntry so a stamped permafrost_zone cannot drift from
   %  the canonical namelist.
   %
   % See also: icemodel.verification.namelists.permafrostzone,
   %  icemodel.verification.setup.validateSurfaceZone

   zone = string(permafrost_zone);
   if zone == "" || ismissing(zone)
      return
   end

   allowed = icemodel.verification.namelists.permafrostzone();
   if ~ismember(zone, allowed)
      error('icemodel:verification:setup:invalidPermafrostZone', ...
         'permafrost_zone "%s" is not in the canonical vocabulary (%s)', ...
         zone, strjoin(allowed, ', '))
   end
end
