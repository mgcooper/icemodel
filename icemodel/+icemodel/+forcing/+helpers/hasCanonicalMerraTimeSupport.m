function tf = hasCanonicalMerraTimeSupport(metadata)
   %HASCANONICALMERRATIMESUPPORT True for the complete MERRA time contract.
   %
   %  tf = icemodel.forcing.helpers.hasCanonicalMerraTimeSupport(metadata)
   %
   % Validates the shared reader/application contract used by MERRA builders,
   % repair tooling, and artifact QA: native coordinates at the reader boundary,
   % center-to-start relabeling only for averaged collections, zero-order hold,
   % and declared 1/1/1/3-hour collection support.

   required = ["merra_source_time_coordinate", ...
      "merra_time_relabel_policy", "merra_time_upsample_policy", ...
      "merra_collection_support_hours"];
   tf = isstruct(metadata) && isscalar(metadata) ...
      && all(isfield(metadata, required));
   if ~tf
      return
   end

   % Reject malformed scalar/string policy values before comparing vocabulary.
   source_coordinate = string(metadata.merra_source_time_coordinate);
   relabel_policy = string(metadata.merra_time_relabel_policy);
   upsample_policy = string(metadata.merra_time_upsample_policy);
   if ~isscalar(source_coordinate) || ~isscalar(relabel_policy) ...
         || ~isscalar(upsample_policy)
      tf = false;
      return
   end

   % Collection support is part of the contract, not optional documentation.
   support = metadata.merra_collection_support_hours;
   support_fields = ["slv", "rad", "flx", "glc"];
   tf = isstruct(support) && isscalar(support) ...
      && all(isfield(support, support_fields));
   if ~tf
      return
   end
   values = zeros(1, numel(support_fields));
   for k = 1:numel(support_fields)
      value = support.(char(support_fields(k)));
      if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
         tf = false;
         return
      end
      values(k) = double(value);
   end
   tf = isequal(values, [1, 1, 1, 3]) ...
      && source_coordinate == "native_at_reader" ...
      && relabel_policy == "time_averaged_center_to_interval_start" ...
      && upsample_policy == "zero_order_hold_over_declared_support";
end
