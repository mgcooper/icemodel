function names = firnCaseManifestFieldNames()
   %FIRNCASEMANIFESTFIELDNAMES Return canonical firn case-manifest fields.
   %
   %  names = icemodel.verification.setup.firnCaseManifestFieldNames()
   %
   % Outputs
   %  names   String array in the exact field order written for each firn
   %          evaluation case.
   %
   % Role
   %  Setup helper used while constructing firn case entries. A firn case is now
   %  a METADATA-ONLY manifest: it records WHICH forcing/eval sources a case
   %  draws on (by source id), the colocation regime, and the comparison
   %  contract - NOT a bundled evaluation.mat/reference.mat data copy. The actual
   %  forcing and eval data live in individual files under the standard icemodel
   %  naming convention (icemodel.forcing.helpers.metfilename / writemet /
   %  writeuserdata): met_<site>_<source>_<window>.mat per forcing source, plus
   %  the per-source Data files for eval. The model runs colocation-agnostically
   %  by picking the met file for the desired source.
   %
   %  Fields
   %    case_id               compact case id ("kanm")
   %    case_type             "firn_observational"
   %    site_id               canonical PROMICE station id ("KAN_M")
   %    site_name             display name
   %    surface_zone          glaciological zone (namelists.surfacezone)
   %    eval_target           capability descriptor array (namelists.evaltarget)
   %    site_location         WGS84 + EPSG:3413 projected coordinates
   %    period                {start, end} the case's evaluation window
   %    forcing_sources       string array of available forcing source ids
   %                          (e.g. ["promice","mar","merra"]) - which met files
   %                          exist to run this site, NOT a bundled data copy
   %    eval_sources          string array of available eval source ids
   %                          (e.g. ["promice_obs","racmo"])
   %    comparison_variables  variables compared at this case
   %    observation_variables observation metadata (e.g. thermistor depths)
   %    colocation            recorded colocation metadata: available models and
   %                          their per-leg windows (METADATA, not bundled data)
   %    native_timestep       observation timestep ("daily")
   %    notes                 provenance note
   %
   % See also: icemodel.verification.setup.caseManifestFieldNames,
   %  icemodel.verification.setup.makeFirnCaseManifestEntry,
   %  icemodel.verification.setup.importPromiceSites

   names = [ ...
      "case_id"
      "case_type"
      "site_id"
      "site_name"
      "surface_zone"
      "eval_target"
      "site_location"
      "period"
      "forcing_sources"
      "eval_sources"
      "comparison_variables"
      "observation_variables"
      "colocation"
      "native_timestep"
      "notes"];
end
