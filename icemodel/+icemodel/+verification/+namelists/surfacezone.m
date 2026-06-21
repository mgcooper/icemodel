function surface_zones = surfacezone()
   %SURFACEZONE Return the supported case-manifest surface-zone values.
   %
   %  surface_zones = icemodel.verification.namelists.surfacezone()
   %
   % Outputs
   %  surface_zones   Supported manifest surface_zone values. These record the
   %                  physical surface regime of a verification case as per-case
   %                  metadata. The regime was formerly encoded by the eval-data
   %                  directory layout (snow/ vs firn/); the taxonomy is now
   %                  dataset-family-flat and the regime lives here instead.
   %
   % Role
   %  Canonical surface-zone vocabulary shared by the setup importers (which
   %  stamp surface_zone onto case manifests) and validators. The single source
   %  of truth for a PROMICE anchor's zone is
   %  icemodel.verification.helpers.promicesiteinfo(site).zone, which uses this
   %  vocabulary. An empty surface_zone ("") is permitted for analytical cases
   %  where the regime is not meaningful (e.g. the Laugh-Tests Colbeck benchmark).

   % Keep this explicit so adding a new zone is a deliberate schema change.
   % Ordered from low to high on a Greenland elevation transect, then the
   % seasonal/dry-snow regimes used by the snow families.
   surface_zones = [ ...
      "bare_ice"
      "ablation"
      "lower_ablation"
      "upper_ablation"
      "percolation"
      "lower_percolation"
      "accumulation"
      "dry_snow"
      "seasonal_snow"];
end
