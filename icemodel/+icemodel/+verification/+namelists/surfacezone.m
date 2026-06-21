function surface_zones = surfacezone()
   %SURFACEZONE Return the supported case-manifest surface-zone values.
   %
   %  surface_zones = icemodel.verification.namelists.surfacezone()
   %
   % Outputs
   %  surface_zones   Supported manifest surface_zone values. These record the
   %                  glaciological substrate REGIME of a verification case as
   %                  per-case metadata - the physical zone the case sits in, not
   %                  the model capability the case exercises (that descriptor is
   %                  eval_target; see icemodel.verification.namelists.evaltarget).
   %
   % Role
   %  Canonical surface-zone vocabulary shared by the setup importers (which
   %  stamp surface_zone onto case manifests) and validators. surface_zone is
   %  the glaciological zone ONLY: where on a Greenland (or off-ice) transect
   %  the case sits. It deliberately does NOT carry capability descriptors such
   %  as "bare_ice" or "seasonal_snow" - those are eval_target values. The single
   %  source of truth for a PROMICE anchor's zone is
   %  icemodel.verification.helpers.promicesiteinfo(site).surface_zone, which uses
   %  this vocabulary. An empty surface_zone ("") is permitted for analytical
   %  cases where the regime is not meaningful (e.g. the Laugh-Tests Colbeck
   %  benchmark).
   %
   % See also: icemodel.verification.namelists.evaltarget,
   %  icemodel.verification.helpers.promicesiteinfo

   % Keep this explicit so adding a new zone is a deliberate schema change.
   % Ordered low-to-high on a Greenland ice-sheet transect (ablation through
   % accumulation), then the off-ice regimes (land, tundra) and "unknown" for
   % non-curated PROMICE stations whose zone is not yet classified.
   surface_zones = [ ...
      "ablation"
      "percolation"
      "wet_snow"
      "dry_snow"
      "accumulation"
      "land"
      "tundra"
      "unknown"];
end
