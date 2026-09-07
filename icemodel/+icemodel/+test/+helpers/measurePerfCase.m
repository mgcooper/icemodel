function perf_data = measurePerfCase(experiment, suite, c, isolation, ...
      config_case, data_root, n_runs, artifact_dir, attempt)
   %MEASUREPERFCASE Measure one case under the selected isolation protocol.
   %
   %  perf_data = icemodel.test.helpers.measurePerfCase(experiment, suite, ...
   %     c, isolation, config_case, data_root, n_runs, artifact_dir, attempt)
   %
   % run_perf_suite and build_perf_baseline both measure through this one
   % function, so the comparison and the accepted baseline always use the
   % same protocol. ISOLATION "process" runs the case in a fresh
   % `matlab -batch` subprocess; "session" times it in this session with
   % `clear functions` hygiene. ATTEMPT names the subprocess spec and
   % result MAT files so retries never overwrite an earlier measurement.
   %
   % See also: icemodel.test.helpers.retryInvalidMeasurement,
   %  icemodel.test.helpers.runPerfCaseSubprocess,
   %  icemodel.test.helpers.runPerfCase

   if isolation == "process"
      % Fresh subprocess per case: no inherited JIT state, persistents,
      % or heap layout. The spec and result MAT files stay in the run's
      % artifact folder so the raw measurement remains auditable.
      spec = struct();
      spec.case_struct = table2struct(c);
      spec.n_runs = n_runs;
      spec.config_case = string(config_case);
      spec.data_root = string(data_root);
      spec.result_file = fullfile(artifact_dir, sprintf( ...
         'isolated_%s_attempt%d_result.mat', c.case_id, attempt));
      spec_file = fullfile(artifact_dir, sprintf( ...
         'isolated_%s_attempt%d_spec.mat', c.case_id, attempt));
      save(spec_file, '-struct', 'spec');

      % The child adds the same source tree this function runs from: this
      % file lives under <source>/+icemodel/+test/+helpers, so four
      % fileparts calls reach the folder that addpath needs.
      source_path = fileparts(fileparts(fileparts(fileparts( ...
         mfilename('fullpath')))));
      matlab_bin = fullfile(matlabroot, 'bin', 'matlab');
      cmd = sprintf(['"%s" -nodisplay -nosplash -batch ' ...
         '"addpath(''%s''); ' ...
         'icemodel.test.helpers.runPerfCaseSubprocess(''%s'')"'], ...
         matlab_bin, source_path, spec_file);
      [status, out] = system(cmd);
      if status ~= 0 || exist(spec.result_file, 'file') ~= 2
         error('icemodel:test:perf:subprocessFailed', ...
            'isolated case %s failed (exit %d):\n%s', ...
            c.case_id, status, out)
      end
      perf_data = load(spec.result_file);
   else
      % In-session hygiene: reset compiled functions and persistents so
      % this case does not inherit the previous case's session state. The
      % experiment's warmup re-warms what the timing needs. The Code
      % Analyzer advisory about this option's performance cost describes
      % the intended effect: uniform cold-ish starts beat fast dirty ones
      % for measurement.
      clear('functions')
      perf_data = icemodel.test.helpers.runPerfCase(experiment, suite, c);
   end
end
