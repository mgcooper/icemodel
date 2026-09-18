function attestation = summarizeMachineState(samples)
   %SUMMARIZEMACHINESTATE Reduce machine-state samples to one attestation.
   %
   %  attestation = icemodel.test.helpers.summarizeMachineState(samples)
   %
   % SAMPLES is the struct array that sampleMachineState produced at run
   % start, after each case, and at run end. The attestation is the record
   % a reader uses to judge a timing run after the fact, and the release
   % snapshot gate reads it through assertReleasePerfBaselineSource.
   %
   % Fields of ATTESTATION
   %  sample_count              Number of samples.
   %  foreign_matlab_processes  Largest foreign MATLAB count seen. NaN when
   %                            any sample could not read the process list.
   %  load_average_at_start     The one-minute load average of the first
   %                            sample, taken before the first measurement.
   %                            The run's own subprocesses raise the later
   %                            samples, so this is the value the release
   %                            gate reads.
   %  load_average_min, load_average_median, load_average_max
   %                            Statistics of the one-minute load average
   %                            over every sample, recorded for the reader.
   %  ac_power                  True only when every sample drew AC power.
   %  load_average_samples, foreign_matlab_samples, ac_power_samples,
   %  sampled_utc               The per-sample values, in order.
   %  probe_errors              Every probe failure text across samples.
   %
   % See also: icemodel.test.helpers.sampleMachineState,
   %  icemodel.test.helpers.assertReleasePerfBaselineSource

   arguments
      samples (:, 1) struct
   end

   loads = reshape([samples.load_average_1min], [], 1);
   foreign = reshape([samples.foreign_matlab_processes], [], 1);
   ac = reshape([samples.ac_power], [], 1);
   errors = arrayfun(@(s) reshape(string(s.probe_errors), [], 1), ...
      samples, 'UniformOutput', false);

   % max and median keep NaN here on purpose: a sample that could not read
   % the machine state must not make the run look quiet.
   attestation = struct( ...
      'sample_count', numel(samples), ...
      'foreign_matlab_processes', max(foreign, [], 'includenan'), ...
      'load_average_at_start', loads(1), ...
      'load_average_min', min(loads, [], 'includenan'), ...
      'load_average_median', median(loads), ...
      'load_average_max', max(loads, [], 'includenan'), ...
      'ac_power', all(ac), ...
      'load_average_samples', loads, ...
      'foreign_matlab_samples', foreign, ...
      'ac_power_samples', ac, ...
      'sampled_utc', reshape([samples.sampled_utc], [], 1), ...
      'probe_errors', vertcat(errors{:}));
end
