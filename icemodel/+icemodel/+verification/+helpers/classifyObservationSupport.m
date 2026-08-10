function support = classifyObservationSupport( ...
      values, field_names, target_field, policy)
   %CLASSIFYOBSERVATIONSUPPORT Apply the flag-support rules to observation rows.
   %
   % support = icemodel.verification.helpers.classifyObservationSupport( ...
   %    values, field_names, target_field, policy)
   %
   % Applies the PROMICE flag rules that decide whether an observation row is
   % supported. Callers pass the observation matrix and the policy; the
   % returned masks are row-shaped logicals they combine as they need.
   %
   % The snow half of the same admission rule is in
   % icemodel.verification.helpers.classifySnowDepth.
   %
   % Inputs
   %  values       - numeric row-by-field observation matrix. Callers holding
   %                 a table extract the columns they need in policy order.
   %  field_names  - string array naming each column of VALUES.
   %  target_field - name of the target observation column. Passed rather than
   %                 read from POLICY because promiceAblationReadiness calls it
   %                 target_field and promiceAblationPolicy, which derives from
   %                 it, exposes the same field as observation_field.
   %  policy       - readiness or ablation policy struct from
   %                 icemodel.verification.namelists.promiceAblationReadiness
   %                 or promiceAblationPolicy. Both carry the same flag lists.
   %
   % Outputs
   %  support - struct of row-shaped logical masks:
   %    target_finite       target observation is finite
   %    quality_finite      every support flag is finite
   %    direct_flags_zero   every direct-zero flag is exactly zero
   %    flag_clean          target_finite & quality_finite & direct_flags_zero
   %    datum_intact        every datum-break flag is finite and zero
   %    gap_flagged         any ordinary-gap flag is finite and nonzero
   %    metadata_flagged    any metadata-only flag is finite and nonzero
   %    station_transition  the station-transition flag is finite and nonzero
   %    unresolved_step     the unresolved-step flag is finite and nonzero
   %
   % Callers compose these masks. The comparator's notion of direct support
   % adds exposed ice, the readiness writer's adds a finite target, and the
   % report builder's adds both.
   %
   % Callers holding a table get the column set from
   % icemodel.verification.helpers.observationSupportFields.
   %
   % See also: icemodel.verification.helpers.observationSupportFields,
   %  icemodel.verification.helpers.classifySnowDepth,
   %  icemodel.verification.namelists.promiceAblationReadiness

   arguments
      values double {mustBeReal}
      field_names string
      target_field (1, 1) string
      policy (1, 1) struct
   end

   % Select every group by name. The policy lists overlap and are ordered
   % differently from each other, so positional indexing would not line each
   % flag up with its meaning.
   support_values = selectColumns(values, field_names, ...
      policy.support_flag_fields);
   direct_values = selectColumns(values, field_names, ...
      policy.direct_zero_flag_fields);
   datum_values = selectColumns(values, field_names, ...
      policy.datum_break_flag_fields);
   gap_values = selectColumns(values, field_names, ...
      policy.ordinary_gap_flag_fields);
   metadata_values = selectColumns(values, field_names, ...
      policy.metadata_only_flag_fields);
   transition_values = selectColumns(values, field_names, ...
      policy.station_transition_flag_field);
   step_values = selectColumns(values, field_names, ...
      policy.unresolved_step_flag_field);
   target_values = selectColumns(values, field_names, target_field);

   support = struct();
   support.target_finite = all(isfinite(target_values), 2);
   support.quality_finite = all(isfinite(support_values), 2);

   % A nonfinite direct-zero flag is not a zero flag, so this test checks
   % finiteness as well. quality_finite covers only the direct-zero flags that
   % are also support flags.
   support.direct_flags_zero = ...
      all(isfinite(direct_values) & direct_values == 0, 2);
   support.flag_clean = support.target_finite ...
      & support.quality_finite & support.direct_flags_zero;
   support.datum_intact = all(isfinite(datum_values) & datum_values == 0, 2);

   % Every set-flag test checks isfinite first: a NaN or an Inf means missing
   % data, not a set flag, and would inflate the counts these masks feed. A
   % finite negative posting is malformed rather than missing and does count as
   % set, which excludes the row.
   support.gap_flagged = anySetFlag(gap_values);
   support.metadata_flagged = anySetFlag(metadata_values);
   support.station_transition = anySetFlag(transition_values);
   support.unresolved_step = anySetFlag(step_values);
end

function selected = selectColumns(values, field_names, wanted)
   %SELECTCOLUMNS Take the named columns, erroring on any that is absent.
   %
   % A missing column would narrow the matrix and make an all() or any()
   % collapse answer a different question than the caller asked.

   wanted = string(wanted);
   [found, where] = ismember(wanted, field_names);
   if ~all(found)
      error('icemodel:verification:classifyObservationSupport:missingField', ...
         'observation matrix has no column named %s', ...
         strjoin(wanted(~found), ', '))
   end
   selected = values(:, where);
end

function flagged = anySetFlag(values)
   %ANYSETFLAG Collapse a flag matrix to one finite nonzero flag per row.
   %
   % Collapsing per row keeps the result a row count when a policy list carries
   % more than one field.

   flagged = any(isfinite(values) & values ~= 0, 2);
end
