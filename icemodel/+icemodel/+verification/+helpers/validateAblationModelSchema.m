function unavailable = validateAblationModelSchema(saved_schema, ...
      current_schema)
   %VALIDATEABLATIONMODELSCHEMA Check a saved cohort still supports the report.
   %
   %  unavailable = ...
   %     icemodel.verification.helpers.validateAblationModelSchema( ...
   %     saved_schema, current_schema)
   %
   % SAVED_SCHEMA is the diagnostic channel list a cohort recorded when it
   % ran, and the list of columns its saved model timetables carry.
   % CURRENT_SCHEMA is the list the running code defines. UNAVAILABLE names
   % the channels the current code defines that the saved cohort does not
   % carry, which the report states rather than treats as an error.
   %
   % Compatibility, not equality. The report reads the channels named by
   % icemodel.verification.namelists.ablationReportChannels. Those must be in
   % the saved results, or the report cannot be built. They must also be in
   % the current namelists. A channel the code does not define has been
   % removed or redefined, so the saved values do not mean what the report
   % says they mean. Every other difference is additive and harmless: a new
   % diagnostic channel does not change a value the cohort already computed.
   %
   % The errors here carry `icemodel:verification:report:` identifiers, not
   % this function's own name. buildAblationEvaluationReport is the only
   % caller, and test_ablation_evaluation_report and test_report_helpers both
   % check for these exact strings. Keep them if this function is renamed.
   %
   % See also: icemodel.verification.namelists.ablationReportChannels,
   %  icemodel.verification.namelists.promiceAblationPolicy,
   %  icemodel.verification.report.buildAblationEvaluationReport

   saved_schema = textSchema(saved_schema, "saved");
   current_schema = textSchema(current_schema, "current");
   consumed = icemodel.verification.namelists.ablationReportChannels();

   % This compares names. A channel that keeps its name while its units or
   % meaning change passes both checks, so a redefined channel is not caught.
   % Bead icemodel-5gs adds the definition stamp that would catch it.
   %
   % Test the current namelists first. A channel they do not define cannot
   % be recovered by running the model again, because a rerun records the same
   % reduced schema. When both lists lack the channel, the rerun advice below
   % would send the owner into a multi-hour loop. This message therefore goes
   % first, and names the two repairs that do work.
   missing_from_current = setdiff(consumed, current_schema, 'stable');
   if ~isempty(missing_from_current)
      error('icemodel:verification:report:incompatibleModelSchema', ...
         ['Report-consumed model channel(s) absent from ' ...
         'icemodel.namelists.budgetoutputs and ' ...
         'icemodel.namelists.cumulativeoutputs: %s. Restore the ' ...
         'channel(s), or drop them from ' ...
         'icemodel.verification.namelists.ablationReportChannels together ' ...
         'with the report code that reads them.'], ...
         strjoin(missing_from_current, ', '))
   end

   % The current code defines the channel, so a cohort that lacks it ran
   % before it existed. A rerun does fix this one.
   missing_from_saved = setdiff(consumed, saved_schema, 'stable');
   if ~isempty(missing_from_saved)
      error('icemodel:verification:report:incompatibleModelSchema', ...
         ['Saved cohort lacks report-consumed model channel(s): %s. ' ...
         'Rerun the cohort.'], strjoin(missing_from_saved, ', '))
   end

   unavailable = setdiff(current_schema, saved_schema, 'stable');
end

function fields = textSchema(fields, role)
   %TEXTSCHEMA Convert one channel list to a text row, or reject it.
   %
   % The strict policy comparison does not cover the channel schema, so a
   % saved MAT file can reach here carrying any type. Reject a value that is
   % not a channel list with invalidAblationPolicy, the identifier every
   % other malformed-policy case raises, not a bare conversion error.

   is_char_cell = iscell(fields) && all(cellfun(@ischar, fields(:)));
   if ~(isstring(fields) || ischar(fields) || is_char_cell)
      error('icemodel:verification:report:invalidAblationPolicy', ...
         'The %s model channel schema must be text', role)
   end
   fields = string(fields(:)');
   if any(ismissing(fields)) || any(strlength(fields) == 0)
      error('icemodel:verification:report:invalidAblationPolicy', ...
         'The %s model channel schema has an empty channel name', role)
   end
end
