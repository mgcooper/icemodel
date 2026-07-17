function [forcing_sources, eval_sources] = colocationSourceLists(colocation)
   %COLOCATIONSOURCELISTS Derive manifest source lists from colocation legs.
   %
   %  [forcing_sources, eval_sources] = ...
   %     icemodel.verification.setup.colocationSourceLists(colocation)
   %
   % Source-list contract
   %  * Observation/protocol legs are eval sources when `eval_staged=true`.
   %    Legacy legs without that field retain their historical `staged`
   %    interpretation:
   %      promice -> promice_obs, sumup -> sumup_obs,
   %      retmip -> retmip_protocol, imau -> imau_obs,
   %      research_site -> research_site_obs, gcnet -> gcnet_obs.
   %  * MAR/MERRA/RACMO staged model legs are eval sources, reported by
   %    product id (for example mar3.11) so multiple versions can coexist.
   %  * Forcing sources must have a staged met file. `forcing_ready` is
   %    advisory runtime metadata: it says whether the met artifact can run
   %    without filling or repair, but it does not suppress source discovery.
   %    RACMO is therefore eval/Data only until a RACMO met source exists.
   %
   % Colocation may also carry metadata-only records such as `anchor`,
   % `source_association`, and `nearest_noncolocated_promice`; those records are
   % provenance, not staged data, and do not create source ids.

   % Observation legs are eval sources. Labels use `_obs` or `_protocol` so
   % they do not collide with native met/model source labels.
   eval_legs = ["promice", "sumup", "retmip", "imau", "research_site", ...
      "gcnet"];
   eval_labels = ["promice_obs", "sumup_obs", "retmip_protocol", ...
      "imau_obs", "research_site_obs", "gcnet_obs"];
   eval_staged = arrayfun(@(source) ...
      evaluationStaged(colocation, source), eval_legs);

   % Model legs are eval sources only when they carry userdata/Data artifacts.
   % Met-only legs can drive forcing but cannot be compared as model data.
   model_legs = icemodel.verification.namelists.rcmsources();
   model_staged = arrayfun(@(source) stagedWithData(colocation, source), ...
      model_legs);
   model_labels = arrayfun(@(source) sourceLabel(colocation, source), ...
      model_legs);
   eval_sources = reshape([eval_labels(eval_staged), ...
      model_labels(model_staged)], [], 1);

   % Forcing sources need a met file. RACMO subsurface remains userdata/Data
   % only, so it is eval-capable but not a forcing source.
   forcing_legs = ["promice", "retmip", "imau", "gcnet", ...
      icemodel.verification.namelists.rcmMetSources()];
   forcing_has_met = arrayfun(@(source) ...
      stagedWithMet(colocation, source), forcing_legs);
   forcing_labels = arrayfun(@(source) sourceLabel(colocation, source), ...
      forcing_legs);
   forcing_sources = reshape(forcing_labels(forcing_has_met), [], 1);
end

function label = sourceLabel(colocation, source)
   %SOURCELABEL Return the manifest source label for one staged leg.
   source = string(source);
   name = char(source);
   if isfield(colocation, name) && isstruct(colocation.(name)) ...
         && isfield(colocation.(name), 'source_id') ...
         && strlength(string(colocation.(name).source_id)) > 0
      label = string(colocation.(name).source_id);
   elseif ismember(source, icemodel.verification.namelists.rcmsources())
      label = icemodel.verification.namelists.rcmProductIds(source);
   else
      label = source;
   end
end

function tf = staged(colocation, source)
   %STAGED True when a colocation leg exists and represents staged data.
   source = char(source);
   tf = false;
   if ~isfield(colocation, source) || ~isstruct(colocation.(source))
      return
   end
   leg = colocation.(source);
   if isfield(leg, 'staged')
      tf = logical(leg.staged);
   end
end

function tf = evaluationStaged(colocation, source)
   %EVALUATIONSTAGED Resolve evaluation availability without conflating runtime.
   source = char(source);
   tf = false;
   if ~isfield(colocation, source) || ~isstruct(colocation.(source))
      return
   end
   leg = colocation.(source);
   if isfield(leg, 'eval_staged')
      % An explicit evaluation fact must override native runtime staging so an
      % RCM-only PROMICE import can expose observations without a false leg.
      tf = logical(leg.eval_staged);
   else
      % Existing manifests predate eval_staged; preserve their source lists.
      tf = staged(colocation, source);
   end
end

function tf = hasMet(leg)
   %HASMET True when a staged leg includes at least one met file.
   tf = isfield(leg, 'met_files') && ~isempty(leg.met_files);
end

function tf = hasData(leg)
   %HASDATA True when a staged leg includes at least one userdata/Data file.
   tf = isfield(leg, 'data_files') && ~isempty(leg.data_files);
end

function tf = stagedWithMet(colocation, source)
   %STAGEDWITHMET True when a staged leg advertises a selectable met artifact.
   source = char(source);
   tf = staged(colocation, source) && hasMet(colocation.(source));
end

function tf = stagedWithData(colocation, source)
   %STAGEDWITHDATA True when a staged leg can be used as an eval source.
   source = char(source);
   tf = staged(colocation, source) && hasData(colocation.(source));
end
