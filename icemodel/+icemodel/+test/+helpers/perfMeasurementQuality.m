function quality = perfMeasurementQuality(case_rows, meta)
   %PERFMEASUREMENTQUALITY Judge a perf run on its measurement quality.
   %
   %  quality = icemodel.test.helpers.perfMeasurementQuality(case_rows, meta)
   %
   % A baseline is a record of measurements, so it is accepted on the
   % quality of those measurements and never on agreement with the baseline
   % it replaces. Trustworthy quality is exactly these five conditions:
   %
   %  samples_valid     every case sample set passed perfSampleValidity;
   %  process_isolation the run measured under isolation="process";
   %  one_machine       the artifact records one nonblank machine identity;
   %  ambient_stable    the ambient anchor was measured and held, and no
   %                    drift override was accepted;
   %  attestation       the machine-state attestation is present.
   %
   % Input
   %  case_rows  Table or struct array with a logical valid column.
   %  meta       Artifact or baseline metadata.
   %
   % Output
   %  quality.passed      True when every condition holds.
   %  quality.conditions  Struct of the five logicals.
   %  quality.reasons     One sentence per failed condition.
   %
   % See also: icemodel.test.helpers.perfSampleValidity,
   %  icemodel.test.helpers.assertReleasePerfBaselineSource,
   %  run_perf_suite, build_perf_baseline

   if istable(case_rows)
      valid = case_rows.valid;
   else
      valid = [case_rows.valid];
   end
   conditions = struct();
   conditions.samples_valid = ~isempty(valid) && all(logical(valid));
   conditions.process_isolation = isfield(meta, 'isolation') ...
      && string(meta.isolation) == "process";
   conditions.one_machine = isfield(meta, 'hostname') ...
      && ~isblanktext(string(meta.hostname));
   drift_accepted = isfield(meta, 'ambient_drift_accepted') ...
      && logical(meta.ambient_drift_accepted);
   % A run that skipped the anchor (run_perf_suite measure_anchor=false)
   % cannot certify ambient stability and says so by name.
   anchor_measured = ~isfield(meta, 'anchor_measured') ...
      || logical(meta.anchor_measured);
   conditions.ambient_stable = isfield(meta, 'ambient_stable') ...
      && logical(meta.ambient_stable) && ~drift_accepted && anchor_measured;
   attestation_fields = ["foreign_matlab_processes", ...
      "load_average_max", "ac_power"];
   conditions.attestation = isfield(meta, 'attestation') ...
      && isstruct(meta.attestation) ...
      && all(isfield(meta.attestation, attestation_fields));

   % One sentence per failed condition, so a reader sees every miss at once.
   messages = struct( ...
      'samples_valid', "at least one case sample set is invalid", ...
      'process_isolation', "the run did not use isolation=""process""", ...
      'one_machine', "the artifact records no machine identity", ...
      'ambient_stable', "the ambient anchor drifted or a drift override was accepted", ...
      'attestation', "the machine-state attestation is missing");
   if ~anchor_measured
      messages.ambient_stable = "the ambient anchor was not measured";
   end
   names = string(fieldnames(conditions));
   failed = names(~cellfun(@(name) conditions.(name), cellstr(names)));
   reasons = arrayfun(@(name) messages.(name), failed);
   quality = struct('passed', isempty(failed), ...
      'conditions', conditions, 'reasons', reshape(reasons, [], 1));
end
