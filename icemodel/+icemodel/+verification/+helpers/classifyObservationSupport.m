function support = classifyObservationSupport( ...
      values, field_names, target_field, policy)
   %CLASSIFYOBSERVATIONSUPPORT Apply the flag-support rules to observation rows.
   %
   % support = icemodel.verification.helpers.classifyObservationSupport( ...
   %    values, field_names, target_field, policy)
   %
   % One owner of how a PROMICE observation row is judged supported. The
   % comparator, the readiness writer, the evaluation runner, and the report
   % builder all consume the same flag lists out of promiceAblationReadiness,
   % and each used to re-derive these masks locally. That let them drift: a
   % second entry in any flag list would have made the readiness writer and the
   % comparator count different row sets for the same site-year, silently,
   % because collapsing a flag matrix to one flag per row is only equivalent to
   % testing it directly while the list has exactly one member.
   %
   % This is the flag counterpart of
   % icemodel.verification.helpers.classifySnowDepth, which already owns the
   % snow half of the same admission rule.
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
   % Callers compose these. The comparator's notion of direct support adds
   % exposed ice, the readiness writer's adds a finite target, and the report
   % builder's adds both; none of those combinations belongs here, because each
   % is a question about a particular consumer rather than about the flags.
   %
   % Callers holding a table get the column set from
   % icemodel.verification.helpers.observationSupportFields rather than
   % rebuilding the union by hand, so adding a flag list to the policy reaches
   % every consumer at once.
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
   % differently from each other, so positional indexing would silently swap
   % two flags' meanings the moment either list changed.
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

   % A nonfinite direct-zero flag is not a zero flag, so the finiteness test
   % has to carry it; quality_finite covers that for the flags that are also
   % support flags, and this keeps the rule true on its own.
   support.direct_flags_zero = ...
      all(isfinite(direct_values) & direct_values == 0, 2);
   support.flag_clean = support.target_finite ...
      & support.quality_finite & support.direct_flags_zero;
   support.datum_intact = all(isfinite(datum_values) & datum_values == 0, 2);

   % isfinite first on every set-flag test. Without it a NaN or an Inf reads
   % as a set flag and inflates the counts these masks feed. A finite negative
   % posting is malformed rather than missing, and does count as set: the
   % conservative direction is to exclude the row rather than admit it.
   support.gap_flagged = anySetFlag(gap_values);
   support.metadata_flagged = anySetFlag(metadata_values);
   support.station_transition = anySetFlag(transition_values);
   support.unresolved_step = anySetFlag(step_values);
end

function selected = selectColumns(values, field_names, wanted)
   %SELECTCOLUMNS Take the named columns, erroring on any that is absent.
   %
   % A missing column would otherwise reduce the matrix silently and make an
   % all() or any() collapse answer a different question than the caller asked.

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
   % Collapsing per row rather than counting the matrix keeps a count a row
   % count as soon as a policy list carries more than one field.

   flagged = any(isfinite(values) & values ~= 0, 2);
end
