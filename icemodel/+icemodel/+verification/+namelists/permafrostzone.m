function permafrost_zones = permafrostzone()
   %PERMAFROSTZONE Return the supported case-manifest permafrost-zone values.
   %
   %  permafrost_zones = icemodel.verification.namelists.permafrostzone()
   %
   % Outputs
   %  permafrost_zones   Supported manifest permafrost_zone values. These record
   %                     the permafrost EXTENT class of the ground a verification
   %                     case sits on, using the standard IPA extent codes
   %                     sampled from the Obu et al. (2019) permafrost-zone map
   %                     (the v1 Brown et al. 1997 source was replaced).
   %
   % Role
   %  Canonical permafrost-zone vocabulary shared by the setup importers (which
   %  stamp permafrost_zone onto case manifests) and validators. permafrost_zone
   %  is ORTHOGONAL to surface_zone: surface_zone records the glaciological
   %  substrate REGIME ("where on the ice/off-ice transect"), permafrost_zone
   %  records the permafrost extent of the GROUND. An ice-sheet case sits on ice,
   %  not permafrost ground, so it carries "none"; an off-ice land/tundra case
   %  carries the Obu extent class sampled at its location. "unknown" is
   %  permitted where a site's permafrost extent has not been resolved.
   %
   %  The non-"none"/"unknown" values are the IPA extent codes:
   %    continuous     90-100% of the area underlain by permafrost (EXTENT "C")
   %    discontinuous  50-90%  (EXTENT "D")
   %    sporadic       10-50%  (EXTENT "S")
   %    isolated       0-10%   (EXTENT "I")
   %
   % See also: icemodel.verification.namelists.surfacezone,
   %  icemodel.verification.helpers.promicesiteinfo

   % Keep this explicit so adding a zone is a deliberate schema change. Ordered
   % most-to-least extensive (continuous through isolated), then "none" for
   % ice-sheet/glacier surfaces that sit on ice rather than permafrost ground,
   % and "unknown" for sites whose permafrost extent is not yet classified.
   permafrost_zones = [ ...
      "continuous"
      "discontinuous"
      "sporadic"
      "isolated"
      "none"
      "unknown"];
end
