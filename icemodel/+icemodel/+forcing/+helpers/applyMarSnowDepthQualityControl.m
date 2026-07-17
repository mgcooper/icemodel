function [T, metadata] = applyMarSnowDepthQualityControl(T, metadata, kwargs)
   %APPLYMARSNOWDEPTHQUALITYCONTROL Mask source-discontinuous SHSN2 years.
   %
   %  [T, metadata] = ... applyMarSnowDepthQualityControl(T)
   %  [T, metadata] = ... applyMarSnowDepthQualityControl(_, metadata)
   %
   % MAR SHSN2 is snow-pack height above ice and is the definition compatible
   % with seasonal snow evaluation. Some selected pixels contain large source
   % resets at annual-file boundaries. This helper screens the daily 00:00
   % samples using an archive-calibrated, scale-relative boundary rule. It
   % masks an interior year only when severe jumps bracket it and its annual
   % median is farther from both neighbours than those neighbours are from
   % each other. One-sided or otherwise ambiguous edges remain finite but are
   % explicitly unverified. SHSN3 is total snow/firn thickness and is never
   % substituted.
   %
   % The default rule is conservative: a boundary jump must exceed both 0.5 m
   % and ten times the adjacent years' pooled 99th-percentile daily increment.
   % The calibration remains explicit so a future source version can be
   % re-audited without changing the provenance vocabulary.
   %
   % See also: icemodel.forcing.buildMarData,
   %  icemodel.verification.auditArtifacts

   arguments
      T timetable
      metadata (1, 1) struct = struct()
      kwargs.jump_ratio (1, 1) double {mustBePositive} = 10
      kwargs.minimum_jump_m (1, 1) double {mustBeNonnegative} = 0.5
   end

   names = string(T.Properties.VariableNames);
   if ~ismember("snowd", names)
      metadata = stampMetadata(metadata, "not_applicable", kwargs, ...
         zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
         zeros(0, 1), 0);
      return
   end

   % A current record is already missing the values needed to recompute the
   % source screen. Re-enforce its durable mask and otherwise remain byte-stable.
   if metadataCurrent(metadata, kwargs)
      masked_years = double(metadata.mar_snowd_masked_years(:));
      remask = ismember(year(T.Time), masked_years) & isfinite(T.snowd);
      if any(remask)
         T.snowd(remask) = NaN;
         metadata.mar_snowd_masked_sample_count = ...
            double(metadata.mar_snowd_masked_sample_count) + nnz(remask);
      end
      return
   end

   % Hour-zero rows are exactly the native daily SHSN2 samples after the
   % support-preserving daily-to-hourly interpolation. Screen only adjacent
   % calendar years with a continuous daily boundary and finite local support.
   daily = hour(T.Time) == 0 & minute(T.Time) == 0 & second(T.Time) == 0;
   years = unique(year(T.Time(daily)))';
   boundary_years = zeros(0, 1);
   jumps = zeros(0, 1);
   references = zeros(0, 1);
   medians = nan(size(years));
   for k = 1:numel(years)
      values = T.snowd(daily & year(T.Time) == years(k));
      medians(k) = median(values, 'omitmissing');
   end
   evaluated = 0;
   for k = 2:numel(years)
      previous_year = years(k - 1);
      current_year = years(k);
      if current_year ~= previous_year + 1
         continue
      end
      previous = T.snowd(daily & year(T.Time) == previous_year);
      current = T.snowd(daily & year(T.Time) == current_year);
      if numel(previous) < 2 || numel(current) < 2 ...
            || ~isfinite(previous(end)) || ~isfinite(current(1))
         continue
      end
      increments = [abs(diff(previous)); abs(diff(current))];
      increments = increments(isfinite(increments));
      if isempty(increments)
         continue
      end
      evaluated = evaluated + 1;
      reference = percentile99(increments);
      jump = abs(current(1) - previous(end));
      if jump > kwargs.minimum_jump_m ...
            && jump > kwargs.jump_ratio * reference
         boundary_years(end + 1, 1) = current_year; %#ok<AGROW>
         jumps(end + 1, 1) = jump; %#ok<AGROW>
         references(end + 1, 1) = reference; %#ok<AGROW>
      end
   end

   % Only a two-sided three-year comparison can identify an isolated source
   % year without assigning blame across a single ambiguous boundary.
   masked_years = zeros(0, 1);
   for k = 2:numel(years) - 1
      has_entry = ismember(years(k), boundary_years);
      has_exit = ismember(years(k + 1), boundary_years);
      if ~has_entry || ~has_exit || any(~isfinite(medians(k - 1:k + 1)))
         continue
      end
      neighbour_distance = abs(medians(k - 1) - medians(k + 1));
      current_distance = min(abs(medians(k) - medians(k - 1)), ...
         abs(medians(k) - medians(k + 1)));
      if current_distance > neighbour_distance
         masked_years(end + 1, 1) = years(k); %#ok<AGROW>
      end
   end
   masked_years = unique(masked_years);

   % Edges adjacent to an isolated masked year are resolved. The year entered
   % by every remaining edge is retained but explicitly unverified, including
   % endpoint cases such as CP1 2019.
   resolved = false(size(boundary_years));
   for k = 1:numel(masked_years)
      resolved = resolved | boundary_years == masked_years(k) ...
         | boundary_years == masked_years(k) + 1;
   end
   unverified_years = boundary_years(~resolved);
   mask = ismember(year(T.Time), masked_years) & isfinite(T.snowd);
   T.snowd(mask) = NaN;
   if evaluated == 0
      status = "insufficient_context";
   else
      status = "applied";
   end
   metadata = stampMetadata(metadata, status, kwargs, boundary_years, ...
      jumps, references, masked_years, unverified_years, nnz(mask));
end

function metadata = stampMetadata(metadata, status, kwargs, ...
      boundary_years, jumps, references, masked_years, unverified_years, ...
      masked_count)
   %STAMPMETADATA Write the stable SHSN2 definition and discontinuity record.
   metadata.mar_snowd_qc_method = ...
      'shsn2_archive_calibrated_year_boundary_screen';
   metadata.mar_snowd_qc_status = char(status);
   metadata.mar_snowd_source = 'SHSN2';
   metadata.mar_snowd_semantics = 'snow_pack_height_above_ice';
   metadata.mar_snowd_shsn3_policy = ...
      'not_used_total_multilayer_snow_firn_thickness';
   metadata.mar_snowd_jump_ratio = kwargs.jump_ratio;
   metadata.mar_snowd_minimum_jump_m = kwargs.minimum_jump_m;
   metadata.mar_snowd_discontinuous_boundary_years = boundary_years(:);
   metadata.mar_snowd_boundary_jump_m = jumps(:);
   metadata.mar_snowd_boundary_reference_p99_m = references(:);
   metadata.mar_snowd_masked_years = masked_years(:);
   metadata.mar_snowd_unverified_years = unverified_years(:);
   metadata.mar_snowd_masked_sample_count = masked_count;
   metadata.mar_snowd_qc_basis = [ ...
      '679 selected-artifact year boundaries showed a ratio gap between ' ...
      '9.87 and 11.8; an isolated interior year is masked only when severe ' ...
      'jumps bracket it and its median is separated from both neighbours; ' ...
      'ambiguous one-sided edges are retained as unverified'];
end

function tf = metadataCurrent(metadata, kwargs)
   %METADATACURRENT Match the exact source definition and calibrated rule.
   required = ["mar_snowd_qc_method", "mar_snowd_jump_ratio", ...
      "mar_snowd_minimum_jump_m", "mar_snowd_masked_years", ...
      "mar_snowd_unverified_years", "mar_snowd_masked_sample_count"];
   tf = all(isfield(metadata, required)) ...
      && string(metadata.mar_snowd_qc_method) ...
      == "shsn2_archive_calibrated_year_boundary_screen" ...
      && double(metadata.mar_snowd_jump_ratio) == kwargs.jump_ratio ...
      && double(metadata.mar_snowd_minimum_jump_m) == kwargs.minimum_jump_m;
end

function value = percentile99(values)
   %PERCENTILE99 Linear-interpolate the finite empirical 99th percentile.
   values = sort(values(:));
   position = 1 + 0.99 * (numel(values) - 1);
   lower = floor(position);
   upper = ceil(position);
   if lower == upper
      value = values(lower);
   else
      weight = position - lower;
      value = values(lower) * (1 - weight) + values(upper) * weight;
   end
end
