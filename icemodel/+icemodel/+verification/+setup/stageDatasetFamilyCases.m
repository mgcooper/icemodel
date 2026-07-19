function [state, alive, skipped] = stageDatasetFamilyCases(items, ...
      empty_state, stage_callback, kwargs)
   %STAGEDATASETFAMILYCASES Stage requested cases with shared skip handling.
   %
   % Importers provide the family-specific staging callback and label callback.
   % This helper owns the common state preallocation, alive mask, skipped-record
   % compaction, skip_missing behavior, and warning shape.

   arguments
      items
      empty_state (1, 1) struct
      stage_callback
      kwargs.skip_missing (1, 1) logical = true
      kwargs.warning_id (1, 1) string
      kwargs.label_callback = []
   end

   n_items = numel(items);
   state = repmat(empty_state, 1, n_items);
   alive = false(1, n_items);
   skipped = repmat(struct('site', "", 'reason', ""), 1, n_items);
   n_skipped = 0;

   for k = 1:n_items
      label = fallbackItemLabel(items, k);
      try
         label = itemLabel(items, k, kwargs.label_callback);
         state(k) = stage_callback(items(k), k);
         alive(k) = true;
      catch err
         if ~kwargs.skip_missing || ~isSkippableMissingDataError(err)
            rethrow(err)
         end
         n_skipped = n_skipped + 1;
         skipped(n_skipped) = struct('site', label, ...
            'reason', string(err.message));
         warning(kwargs.warning_id, 'skipping %s: %s', label, err.message);
      end
   end

   skipped = skipped(1:n_skipped);
end

function label = fallbackItemLabel(items, k)
   %FALLBACKITEMLABEL Keep skipped records useful if label lookup itself fails.
   try
      label = string(items(k));
   catch
      label = "item " + k;
   end
end

function label = itemLabel(items, k, label_callback)
   %ITEMLABEL Resolve the case/site label used in skipped records.
   if ~isempty(label_callback)
      label = string(label_callback(items(k), k));
   else
      label = string(items(k));
   end
end

function tf = isSkippableMissingDataError(err)
   %ISSKIPPABLEMISSINGDATAERROR True for absent source/cache/window failures.
   id = string(err.identifier);
   skippable = [ ...
      "icemodel:verification:importRetmip:missingSurfaceFile"
      "icemodel:verification:importRetmip:missingNativeSource"
      "icemodel:verification:importRetmip:missingGcnetVandecrux"
      "icemodel:verification:importImau:missingHourlySource"
      "icemodel:verification:importImau:missingDailyQa"
      "icemodel:verification:importLaughTests:missingSource"
      "icemodel:verification:importSumup:missingSource"
      "icemodel:verification:importSumup:missingObservation"
      "icemodel:verification:importResearchSites:missingObservation"
      "icemodel:verification:fetchSumup:missingSources"
      "icemodel:forcing:readPromiceAws:sourceNotFound"
      "icemodel:forcing:readPromiceAws:stationNotFound"];
   tf = any(id == skippable) ...
      || endsWith(id, ":fileNotFound") ...
      || endsWith(id, ":emptyWindow") ...
      || (contains(id, ":empty") && endsWith(id, "Window"));
end
