function validateSurfaceZone(surface_zone)
   %VALIDATESURFACEZONE Validate a case-manifest surface_zone value.
   %
   %  icemodel.verification.setup.validateSurfaceZone(surface_zone)
   %
   % Inputs
   %  surface_zone   The surface_zone value from a case-manifest entry. An
   %                 empty value ("") is permitted (analytical cases where the
   %                 surface regime is not meaningful); any non-empty value
   %                 must be a member of the canonical vocabulary published by
   %                 icemodel.verification.namelists.surfacezone.
   %
   % Role
   %  Setup-side schema gate shared by makeCaseManifestEntry and
   %  makeFirnCaseManifestEntry so a stamped surface_zone cannot drift from the
   %  canonical namelist.
   %
   % See also: icemodel.verification.namelists.surfacezone

   zone = string(surface_zone);
   if zone == "" || ismissing(zone)
      return
   end

   allowed = icemodel.verification.namelists.surfacezone();
   if ~ismember(zone, allowed)
      error('icemodel:verification:setup:invalidSurfaceZone', ...
         'surface_zone "%s" is not in the canonical vocabulary (%s)', ...
         zone, strjoin(allowed, ', '))
   end
end
