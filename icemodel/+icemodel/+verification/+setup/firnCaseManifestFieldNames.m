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
   %  Setup helper used while constructing firn case entries. A firn case has a
   %  FORCING-AGNOSTIC manifest: the eval target IS bundled as a data-only
   %  observations.mat referenced via evaluation_file (same contract as
   %  ESM-SnowMIP/SUMup), but the FORCING is not bundled or stipulated. The
   %  manifest records WHICH forcing/eval sources a case draws on (by source id,
   %  INFORMATIONAL), the colocation regime, and the comparison contract. The
   %  forcing data lives in individual files under the standard icemodel naming
   %  convention (icemodel.forcing.helpers.metfilename / writemet / writeuserdata):
   %  met_<site>_<source>_<window>.mat per forcing source, plus the per-source
   %  Data files. The model runs colocation-agnostically by picking the met file
   %  for the desired source - any forcing usable at runtime without rewriting
   %  observations.mat.
   %
   %  Fields
   %    case_id               compact case id ("kanm")
   %    case_type             "firn_observational"
   %    site_id               canonical PROMICE station id ("KAN_M")
   %    site_name             display name
   %    surface_zone          glaciological zone (namelists.surfacezone)
   %    eval_target           capability descriptor array (namelists.evaltarget)
   %    permafrost_zone       permafrost extent class, ORTHOGONAL to surface_zone
   %                          (namelists.permafrostzone): off-ice land/tundra
   %                          sites carry the Obu et al. (2019) extent, ice-sheet
   %                          /glacier sites carry "none".
   %    site_location         WGS84 + EPSG:3413 projected coordinates
   %    period                {start, end} the case's evaluation window
   %    evaluation_file       case-relative path to the data-only observations.mat
   %                          (the eval target bundle), e.g. "kanm/observations.mat".
   %                          This is the forcing-AGNOSTIC eval artifact comparecase
   %                          loads; forcing is discovered at runtime, not here.
   %    forcing_sources       string array of forcing source ids staged alongside
   %                          (INFORMATIONAL only - NOT load-bearing; any forcing
   %                          file may be used at runtime regardless of this list)
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
      "permafrost_zone"
      "site_location"
      "period"
      "evaluation_file"
      "forcing_sources"
      "eval_sources"
      "comparison_variables"
      "observation_variables"
      "colocation"
      "native_timestep"
      "notes"];
end
