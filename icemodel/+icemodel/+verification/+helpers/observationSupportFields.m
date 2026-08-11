function fields = observationSupportFields(target_field, policy)
   %OBSERVATIONSUPPORTFIELDS Columns the observation support rules read.
   %
   % fields = icemodel.verification.helpers.observationSupportFields( ...
   %    target_field, policy)
   %
   % Returns the columns
   % icemodel.verification.helpers.classifyObservationSupport reads.
   %
   % promiceAblationPolicy derives required_observation_fields from this, so
   % the comparator gets the set through the policy; the readiness writer, the
   % evaluation runner, and the report builder call it directly.
   %
   % The target column comes first, then every flag group in policy order.
   % Duplicates are dropped because the policy lists overlap.
   %
   % Inputs
   %  target_field - name of the target observation column.
   %                 promiceAblationReadiness calls this target_field and
   %                 promiceAblationPolicy, which derives from it, exposes the
   %                 same field as observation_field, so callers pass their own.
   %  policy       - readiness or ablation policy struct.
   %
   % Outputs
   %  fields - string array of column names.
   %
   % See also: icemodel.verification.helpers.classifyObservationSupport

   arguments
      target_field (1, 1) string
      policy (1, 1) struct
   end

   fields = unique([target_field, ...
      string(policy.support_flag_fields), ...
      string(policy.direct_zero_flag_fields), ...
      string(policy.datum_break_flag_fields), ...
      string(policy.ordinary_gap_flag_fields), ...
      string(policy.metadata_only_flag_fields), ...
      string(policy.station_transition_flag_field), ...
      string(policy.unresolved_step_flag_field)], 'stable');
end
