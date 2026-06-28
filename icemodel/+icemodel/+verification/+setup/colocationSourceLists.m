function [forcing_sources, eval_sources] = colocationSourceLists(colocation)
   %COLOCATIONSOURCELISTS Derive source ids from staged colocation legs.
   %
   %  [forcing_sources, eval_sources] = ...
   %     icemodel.verification.setup.colocationSourceLists(colocation)
   %
   % Role
   %  Shared manifest helper for source-list metadata. `forcing_sources` names
   %  staged sources with valid met forcing files; `eval_sources` names every
   %  staged observation/model source that can be used as a comparison target.

   % Observation legs are eval sources. Labels use `_obs` or `_protocol` so
   % they do not collide with native met/model source labels.
   eval_legs = ["promice", "sumup", "retmip", "imau", "research_site", ...
      "gcnet"];
   eval_labels = ["promice_obs", "sumup_obs", "retmip_protocol", ...
      "imau_obs", "research_site_obs", "gcnet_obs"];
   eval_staged = arrayfun(@(source) staged(colocation, source), eval_legs);

   % Model legs are eval sources even when they are also forcing sources. This
   % supports model-model comparisons, for example MERRA-forced vs staged MAR.
   model_legs = ["mar", "merra", "racmo"];
   model_staged = arrayfun(@(source) staged(colocation, source), model_legs);
   eval_sources = reshape([eval_labels(eval_staged), ...
      model_legs(model_staged)], [], 1);

   % Forcing sources need a met file. RACMO subsurface remains userdata/Data
   % only, so it is eval-capable but not a forcing source.
   forcing_legs = ["promice", "retmip", "imau", "gcnet", "mar", "merra"];
   forcing_has_met = arrayfun(@(source) ...
      stagedWithMet(colocation, source), forcing_legs);
   forcing_sources = reshape(forcing_legs(forcing_has_met), [], 1);
end

function tf = staged(colocation, source)
   %STAGED True when a colocation leg exists and is marked staged.
   source = char(source);
   tf = isfield(colocation, source) ...
      && isstruct(colocation.(source)) ...
      && isfield(colocation.(source), 'staged') ...
      && logical(colocation.(source).staged);
end

function tf = hasMet(leg)
   %HASMET True when a staged leg includes at least one met file.
   tf = isfield(leg, 'met_files') && ~isempty(leg.met_files);
end

function tf = stagedWithMet(colocation, source)
   %STAGEDWITHMET True when a staged leg can be used as a forcing source.
   source = char(source);
   tf = staged(colocation, source) && hasMet(colocation.(source));
end
