function assertReleasePerfBaselineSource(PerfBaseline, meta, kwargs)
   %ASSERTRELEASEPERFBASELINESOURCE Verify a rolling perf source for release.
   %
   %  icemodel.test.helpers.assertReleasePerfBaselineSource(PerfBaseline, meta)
   %
   % snapshot_perf_baseline copies the rolling file into a frozen release
   % file and does not re-measure, so this check reads the recorded
   % metadata instead. A release source must pass the five quality
   % conditions of perfMeasurementQuality, must have been measured on the
   % machine that runs the snapshot, and must carry an attestation with zero
   % foreign MATLAB processes, a maximum one-minute load average at or below
   % perfMeasurementPolicy().load_average_max, and AC power. A rolling
   % baseline that misses any of these may exist; it cannot become a release
   % baseline.
   %
   % Name-value
   %  current_identity  Identity of the machine running the snapshot.
   %                    Defaults to machineHostname(). Tests inject it.
   %  policy            Threshold source. Defaults to perfMeasurementPolicy().
   %
   % See also: icemodel.test.helpers.perfMeasurementQuality,
   %  icemodel.test.helpers.snapshotBaseline, snapshot_perf_baseline

   arguments
      PerfBaseline
      meta (1, 1) struct
      kwargs.current_identity (1, 1) string = ...
         icemodel.test.helpers.machineHostname()
      kwargs.policy (1, 1) struct = ...
         icemodel.test.helpers.perfMeasurementPolicy()
   end

   quality = icemodel.test.helpers.perfMeasurementQuality(PerfBaseline, meta);
   if ~quality.passed
      error('icemodel:test:releasePerfSourceQuality', ...
         ['The rolling perf baseline fails the measurement quality ', ...
         'conditions: %s.'], char(strjoin(quality.reasons, "; ")))
   end

   % Internal consistency is not enough: a file measured entirely on another
   % machine is self-consistent, so the identity must equal this machine's.
   source_identity = icemodel.test.helpers.normalizeMachineIdentity( ...
      string(meta.hostname));
   if source_identity ~= kwargs.current_identity
      error('icemodel:test:releasePerfSourceForeignHost', ...
         ['The rolling perf baseline was measured on "%s"; this ', ...
         'machine is "%s".'], char(source_identity), ...
         char(kwargs.current_identity))
   end

   attestation = meta.attestation;
   if ~(attestation.foreign_matlab_processes == 0)
      error('icemodel:test:releasePerfSourceForeignMatlab', ...
         ['The rolling perf baseline saw %g foreign MATLAB processes ', ...
         'during measurement; a release source needs zero.'], ...
         attestation.foreign_matlab_processes)
   end
   if ~(attestation.load_average_max <= kwargs.policy.load_average_max)
      error('icemodel:test:releasePerfSourceLoadAverage', ...
         ['The rolling perf baseline saw a one-minute load average of ', ...
         '%.2f; a release source needs at most %.2f.'], ...
         attestation.load_average_max, kwargs.policy.load_average_max)
   end
   if ~attestation.ac_power
      error('icemodel:test:releasePerfSourceBatteryPower', ...
         ['The rolling perf baseline was measured on battery power; a ', ...
         'release source needs AC power throughout.'])
   end
end
