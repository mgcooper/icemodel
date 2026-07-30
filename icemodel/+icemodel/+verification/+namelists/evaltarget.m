function eval_targets = evaltarget()
   %EVALTARGET Return the supported case-manifest eval-target descriptors.
   %
   %  eval_targets = icemodel.verification.namelists.evaltarget()
   %
   % Outputs
   %  eval_targets   Supported manifest eval_target values. These name which
   %                 model CAPABILITY a verification case exercises - the
   %                 process the case is staged to evaluate - independent of the
   %                 glaciological zone the case sits in (that is surface_zone).
   %
   % Role
   %  Canonical eval-target vocabulary shared by the setup importers (which stamp
   %  eval_target onto case manifests) and validators. eval_target is a string
   %  ARRAY: one case can exercise several capabilities (e.g. a KAN ablation site
   %  that carries both a seasonal snowpack and a bare-ice melt season exercises
   %  both "seasonal_snow" and "bare_ice"). It answers "what does this case test?"
   %  whereas surface_zone answers "where is this case?". Examples:
   %    KAN_L surface_zone="ablation",    eval_target=["seasonal_snow","bare_ice"]
   %    KAN_U surface_zone="percolation", eval_target=["seasonal_snow","firn"]
   %  An empty eval_target (string(0,1)) is permitted for cases that do not
   %  exercise a curated capability (e.g. analytical Colbeck benchmarks).
   %
   % See also: icemodel.verification.namelists.surfacezone,
   %  icemodel.verification.setup.promiceSiteCatalog

   % Keep this explicit so adding a capability descriptor is a deliberate schema
   % change. "seasonal_snow" / "bare_ice" / "firn" / "ablation" are the capability
   % axes the firn and snow lanes exercise; "land"/"tundra" surfaces use no
   % glaciological capability and carry an empty eval_target.
   eval_targets = [ ...
      "seasonal_snow"
      "bare_ice"
      "firn"
      "ablation"];
end
