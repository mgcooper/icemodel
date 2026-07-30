function [Data, metadata] = applyRacmoPrecipitationQualityControl(Data, prior_metadata)
   %APPLYRACMOPRECIPITATIONQUALITYCONTROL Enforce nonnegative RACMO ppt.
   %
   %  [Data, metadata] = ... applyRacmoPrecipitationQualityControl(Data)
   %  [Data, metadata] = ... applyRacmoPrecipitationQualityControl( ...
   %     Data, prior_metadata)
   %
   % RACMO labels `precip` as a precipitation flux, which is physically
   % nonnegative, but the native model field contains small negative numerical
   % undershoots. Apply the invariant after spatial sampling/remapping and
   % temporal interpolation, when the builder has converted the public `ppt`
   % channel to canonical m s-1. Every finite negative sample becomes exactly
   % zero; missing values, legitimate positives, the time axis, and unrelated
   % channels are preserved.
   %
   % The returned flat metadata fields are suitable for both freshly built
   % artifacts and the exact-reference repair path. There is deliberately no
   % magnitude threshold: precipitation cannot be negative, while QA retains
   % the input minimum and replacement count as provenance. Pass the prior flat
   % QC contract when repairing an already-canonical artifact so the second pass
   % preserves the original input minimum/count as well as the data.
   %
   % See also: icemodel.forcing.buildRacmoData

   arguments
      Data timetable
      prior_metadata struct = struct()
   end

   % Define the stable metadata contract even when a reduced source omits ppt.
   metadata = struct( ...
      'racmo_ppt_qc_method', "negative_to_zero", ...
      'racmo_ppt_qc_stage', ...
      "after_spatial_sampling_and_temporal_interpolation", ...
      'racmo_ppt_qc_source_variable', "precip", ...
      'racmo_ppt_qc_basis', ...
      "RACMO precip is a physical precipitation flux; finite negative " + ...
      "model/sampling/interpolation undershoot is zero", ...
      'racmo_ppt_qc_status', "not_applicable", ...
      'racmo_ppt_qc_replaced_count', 0, ...
      'racmo_ppt_qc_input_minimum', NaN, ...
      'racmo_ppt_qc_output_minimum', NaN);

   names = string(Data.Properties.VariableNames);
   if ~ismember("ppt", names)
      return
   end

   % Record the source-finalized input distribution before enforcing the
   % physical invariant, preserving all nonfinite values without fabrication.
   values = Data.ppt;
   finite = isfinite(values);
   if any(finite, 'all')
      metadata.racmo_ppt_qc_input_minimum = min(values(finite), [], 'all');
   end
   negative = finite & values < 0;
   values(negative) = 0;
   Data.ppt = values;

   % Report whether this pass repaired anything. Reapplying the helper is a
   % data-level no-op; passing the prior contract also preserves metadata.
   metadata.racmo_ppt_qc_replaced_count = nnz(negative);
   if any(negative, 'all')
      metadata.racmo_ppt_qc_status = "applied";
   else
      metadata.racmo_ppt_qc_status = "unchanged";
   end
   finite = isfinite(values);
   if any(finite, 'all')
      metadata.racmo_ppt_qc_output_minimum = min(values(finite), [], 'all');
   end

   % A second repair must retain the original input minimum and replacement
   % count, not rewrite provenance as if the source never contained negatives.
   % Preserve a complete compatible prior contract only when this pass makes no
   % new replacements and the already-repaired output remains nonnegative.
   if ~any(negative, 'all') && compatiblePriorMetadata(prior_metadata, metadata)
      fields = fieldnames(metadata);
      for k = 1:numel(fields)
         metadata.(fields{k}) = prior_metadata.(fields{k});
      end
   end
end

function tf = compatiblePriorMetadata(prior, current)
   %COMPATIBLEPRIORMETADATA Validate preserved idempotent QC provenance.
   fields = fieldnames(current);
   tf = isscalar(prior) && all(isfield(prior, fields));
   if ~tf
      return
   end

   % Reject non-scalar or non-text contract tokens before converting them to
   % strings, so malformed prior metadata falls back cleanly instead of erroring.
   text_fields = { ...
      'racmo_ppt_qc_method', 'racmo_ppt_qc_stage', ...
      'racmo_ppt_qc_source_variable', 'racmo_ppt_qc_basis', ...
      'racmo_ppt_qc_status'};
   for k = 1:numel(text_fields)
      value = prior.(text_fields{k});
      valid_text = (ischar(value) && isrow(value)) ...
         || (isstring(value) && isscalar(value));
      if ~valid_text
         tf = false;
         return
      end
   end

   % Require scalar, internally consistent status/count/range provenance before
   % comparing the settled method, stage, source, and physical basis.
   status = string(prior.racmo_ppt_qc_status);
   count = prior.racmo_ppt_qc_replaced_count;
   input_minimum = prior.racmo_ppt_qc_input_minimum;
   output_minimum = prior.racmo_ppt_qc_output_minimum;
   tf = isscalar(status) && ismember(status, ["applied", "unchanged"]) ...
      && isnumeric(count) && isscalar(count) && isfinite(count) ...
      && count >= 0 && count == fix(count) ...
      && isnumeric(input_minimum) && isscalar(input_minimum) ...
      && isfinite(input_minimum) ...
      && isnumeric(output_minimum) && isscalar(output_minimum) ...
      && isfinite(output_minimum) && output_minimum >= 0;
   if ~tf
      return
   end

   % Applied provenance must describe at least one negative replacement;
   % unchanged provenance must describe an already nonnegative input.
   if status == "applied"
      tf = count > 0 && input_minimum < 0;
   else
      tf = count == 0 && input_minimum >= 0;
   end
   if ~tf
      return
   end

   % Preserve provenance only for this exact contract and current output
   % minimum. A different physical basis or changed payload gets fresh metadata.
   tf = string(prior.racmo_ppt_qc_method) ...
      == string(current.racmo_ppt_qc_method) ...
      && string(prior.racmo_ppt_qc_stage) ...
      == string(current.racmo_ppt_qc_stage) ...
      && string(prior.racmo_ppt_qc_source_variable) ...
      == string(current.racmo_ppt_qc_source_variable) ...
      && string(prior.racmo_ppt_qc_basis) ...
      == string(current.racmo_ppt_qc_basis) ...
      && output_minimum == current.racmo_ppt_qc_output_minimum;
end
