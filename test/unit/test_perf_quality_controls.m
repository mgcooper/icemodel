function tests = test_perf_quality_controls
   %TEST_PERF_QUALITY_CONTROLS Test the perf measurement quality controls.
   %
   % Covers the machine-state attestation (sampleMachineState and
   % summarizeMachineState), the five-condition quality verdict
   % (perfMeasurementQuality), and the release-only source check that
   % snapshot_perf_baseline applies through snapshotBaseline
   % (assertReleasePerfBaselineSource). Every probe is injected, so no test
   % reads the live machine or runs a model.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   %SETUPONCE Put the snapshot tools on the path for the refusal test.
   root = icemodel.internal.fullpath();
   original_path = path;
   testCase.addTeardown(@() path(original_path));
   addpath(fullfile(root, 'test', 'tools'));
end

function runner = fakeProbeRunner(ps_text, uptime_text, pmset_text)
   %FAKEPROBERUNNER Return a command runner that answers the three probes.
   runner = @(command) answer(command);
   function [status, text] = answer(command)
      status = 0;
      if startsWith(command, "ps")
         text = char(ps_text);
      elseif startsWith(command, "uptime")
         text = char(uptime_text);
      else
         text = char(pmset_text);
      end
   end
end

function text = fakePsText(rows)
   %FAKEPSTEXT Render a ps -axo pid,ppid,command listing.
   lines = ["  PID  PPID COMMAND"; rows(:)];
   text = strjoin(lines, newline);
end

function test_sample_counts_only_foreign_matlab_processes(testCase)
   % A MATLAB process that descends from this run is not foreign; every
   % other MATLAB binary is, and MCP servers that only name MATLAB in their
   % path are not MATLAB.
   binary = "/Applications/MATLAB_R2025b.app/bin/maca64/MATLAB";
   rows = [ ...
      "  100     1 " + binary + " -nodesktop"; ...              % own
      "  200   100 /bin/sh matlab"; ...                          % own child
      "  201   200 " + binary + " -batch case"; ...              % own grandchild
      "  300     1 " + binary + " -desktop"; ...                 % foreign
      "  301   300 " + binary + " -batch other"; ...             % foreign child
      "  400     1 /Users/x/mcp/bin/matlab-mcp-core-server --matlab-root /Applications/MATLAB_R2025b.app"];
   runner = fakeProbeRunner(fakePsText(rows), ...
      "22:26  up 60 days, load averages: 4.24 3.77 3.24", ...
      "Now drawing from 'AC Power'");
   returned = icemodel.test.helpers.sampleMachineState( ...
      command_runner=runner, own_pid=100, platform="mac");
   testCase.verifyEqual(returned.foreign_matlab_processes, 2);
   testCase.verifyEqual(returned.load_average_1min, 4.24);
   testCase.verifyTrue(returned.ac_power);
   testCase.verifyEmpty(returned.probe_errors);
   testCase.verifyEqual(returned.sampled_utc.TimeZone, 'UTC');
end

function test_sample_parses_linux_uptime_and_battery_power(testCase)
   % The Linux uptime format has a comma list, and battery power reads false
   % from pmset on macOS.
   runner = fakeProbeRunner(fakePsText(strings(0, 1)), ...
      " 10:00:00 up 3 days, load average: 0.52, 0.60, 0.70", ...
      "Now drawing from 'Battery Power'");
   returned = icemodel.test.helpers.sampleMachineState( ...
      command_runner=runner, own_pid=1, platform="mac");
   testCase.verifyEqual(returned.foreign_matlab_processes, 0);
   testCase.verifyEqual(returned.load_average_1min, 0.52);
   testCase.verifyFalse(returned.ac_power);
end

function test_sample_records_probe_failures_as_nan(testCase)
   % A failed probe must not make the machine look quiet: the field is NaN
   % and the failure text is kept.
   runner = @(~) deal(1, 'permission denied');
   returned = icemodel.test.helpers.sampleMachineState( ...
      command_runner=runner, own_pid=1, platform="mac");
   testCase.verifyTrue(isnan(returned.foreign_matlab_processes));
   testCase.verifyTrue(isnan(returned.load_average_1min));
   testCase.verifyFalse(returned.ac_power);
   testCase.verifyEqual(numel(returned.probe_errors), 3);
end

function test_sample_reads_linux_power_supply_tree(testCase)
   % Off macOS, a machine without a battery entry is on mains power, and a
   % machine with one reads the adapter's online flag.
   fixture = testCase.applyFixture( ...
      matlab.unittest.fixtures.TemporaryFolderFixture);
   runner = fakeProbeRunner(fakePsText(strings(0, 1)), ...
      "load average: 1.00, 1.00, 1.00", "");
   returned = icemodel.test.helpers.sampleMachineState( ...
      command_runner=runner, own_pid=1, platform="linux", ...
      power_supply_dir=fixture.Folder);
   testCase.verifyTrue(returned.ac_power);

   mkdir(fullfile(fixture.Folder, 'BAT0'));
   mkdir(fullfile(fixture.Folder, 'AC0'));
   writelines("0", fullfile(fixture.Folder, 'AC0', 'online'));
   returned = icemodel.test.helpers.sampleMachineState( ...
      command_runner=runner, own_pid=1, platform="linux", ...
      power_supply_dir=fixture.Folder);
   testCase.verifyFalse(returned.ac_power);
   writelines("1", fullfile(fixture.Folder, 'AC0', 'online'));
   returned = icemodel.test.helpers.sampleMachineState( ...
      command_runner=runner, own_pid=1, platform="linux", ...
      power_supply_dir=fixture.Folder);
   testCase.verifyTrue(returned.ac_power);
end

function runner = fakeWindowsRunner(process_text, queue_text, battery_text)
   %FAKEWINDOWSRUNNER Return a runner that answers the three PowerShell probes.
   runner = @(command) answer(command);
   function [status, text] = answer(command)
      status = 0;
      if contains(command, "Win32_Process")
         text = char(process_text);
      elseif contains(command, "ProcessorQueueLength")
         text = char(queue_text);
      else
         text = char(battery_text);
      end
   end
end

function test_sample_reads_windows_probes(testCase)
   % Windows has no ps, uptime, or pmset: the process list comes from
   % Win32_Process in the same "pid ppid command" shape, the load stand-in
   % is the processor queue length, and Win32_Battery gives the power state.
   binary = "C:\Program Files\MATLAB\R2025b\bin\win64\MATLAB.exe";
   lowercase_binary = "c:\program files\matlab\r2025b\bin\win64\matlab.exe";
   rows = [ ...
      "100 1 """ + binary + """ -nodesktop"; ...
      "150 100"; ...                                  % own child, no command
      "201 150 """ + binary + """ -batch case"; ...    % own grandchild
      "300 1 """ + binary + """"; ...                  % foreign
      "301 1 """ + lowercase_binary + """ -batch x"; ... % foreign, lowercase
      "400 1 C:\Windows\explorer.exe"];
   runner = fakeWindowsRunner(strjoin(rows, newline), "3", "2");
   returned = icemodel.test.helpers.sampleMachineState( ...
      command_runner=runner, own_pid=100, platform="windows");
   testCase.verifyEqual(returned.foreign_matlab_processes, 2);
   testCase.verifyEqual(returned.load_average_1min, 3);
   testCase.verifyTrue(returned.ac_power);
   testCase.verifyEmpty(returned.probe_errors);

   % A desktop without a battery draws mains power; two batteries must both
   % report an AC-connected state; a discharging battery does not; a blank
   % queue reading is an unknown load.
   runner = fakeWindowsRunner(strjoin(rows, newline), "0", "none");
   returned = icemodel.test.helpers.sampleMachineState( ...
      command_runner=runner, own_pid=100, platform="windows");
   testCase.verifyTrue(returned.ac_power);
   runner = fakeWindowsRunner(strjoin(rows, newline), "0", "2" + newline + "3");
   returned = icemodel.test.helpers.sampleMachineState( ...
      command_runner=runner, own_pid=100, platform="windows");
   testCase.verifyTrue(returned.ac_power);
   runner = fakeWindowsRunner(strjoin(rows, newline), "0", "2" + newline + "1");
   returned = icemodel.test.helpers.sampleMachineState( ...
      command_runner=runner, own_pid=100, platform="windows");
   testCase.verifyFalse(returned.ac_power);
   runner = fakeWindowsRunner(strjoin(rows, newline), "", "1");
   returned = icemodel.test.helpers.sampleMachineState( ...
      command_runner=runner, own_pid=100, platform="windows");
   testCase.verifyFalse(returned.ac_power);
   testCase.verifyTrue(isnan(returned.load_average_1min));
   testCase.verifyEqual(numel(returned.probe_errors), 1);
end

function samples = fakeSamples(loads, foreign, ac)
   %FAKESAMPLES Build a sample array without probing the machine.
   n = numel(loads);
   samples = repmat(struct('foreign_matlab_processes', 0, ...
      'load_average_1min', 0, 'ac_power', true, ...
      'sampled_utc', datetime('now', 'TimeZone', 'UTC'), ...
      'probe_errors', strings(0, 1)), n, 1);
   for k = 1:n
      samples(k).foreign_matlab_processes = foreign(k);
      samples(k).load_average_1min = loads(k);
      samples(k).ac_power = ac(k);
   end
end

function test_summary_reduces_samples_to_the_attestation(testCase)
   % The attestation keeps the worst foreign count, the load statistics,
   % and AC power only when every sample drew it.
   samples = fakeSamples([1.0 3.5 2.0], [0 1 0], [true true true]);
   returned = icemodel.test.helpers.summarizeMachineState(samples);
   testCase.verifyEqual(returned.sample_count, 3);
   testCase.verifyEqual(returned.foreign_matlab_processes, 1);
   testCase.verifyEqual(returned.load_average_at_start, 1.0);
   testCase.verifyEqual(returned.load_average_min, 1.0);
   testCase.verifyEqual(returned.load_average_median, 2.0);
   testCase.verifyEqual(returned.load_average_max, 3.5);
   testCase.verifyTrue(returned.ac_power);
   testCase.verifyEqual(returned.load_average_samples, [1.0; 3.5; 2.0]);

   samples = fakeSamples([1.0 NaN], [0 NaN], [true false]);
   returned = icemodel.test.helpers.summarizeMachineState(samples);
   testCase.verifyTrue(isnan(returned.foreign_matlab_processes));
   testCase.verifyTrue(isnan(returned.load_average_max));
   testCase.verifyFalse(returned.ac_power);
end

function [rows, meta] = compliantSource()
   %COMPLIANTSOURCE Build rows and metadata that pass every condition.
   rows = table([true; true], [10; 20], 'VariableNames', ...
      {'valid', 'median_wall_s'});
   meta = struct('isolation', "process", 'hostname', "macbook-air-2", ...
      'ambient_stable', true, 'ambient_drift_accepted', false, ...
      'attestation', icemodel.test.helpers.summarizeMachineState( ...
      fakeSamples([1.5 2.0 1.0], [0 0 0], [true true true])));
end

function test_quality_passes_only_on_all_five_conditions(testCase)
   % Each condition fails on its own, and the comparison verdict is not
   % one of them.
   [rows, meta] = compliantSource();
   returned = icemodel.test.helpers.perfMeasurementQuality(rows, meta);
   testCase.verifyTrue(returned.passed);
   testCase.verifyEmpty(returned.reasons);

   invalid_rows = rows;
   invalid_rows.valid(2) = false;
   returned = icemodel.test.helpers.perfMeasurementQuality(invalid_rows, meta);
   testCase.verifyFalse(returned.passed);
   testCase.verifyFalse(returned.conditions.samples_valid);

   session_meta = meta;
   session_meta.isolation = "session";
   returned = icemodel.test.helpers.perfMeasurementQuality(rows, session_meta);
   testCase.verifyFalse(returned.conditions.process_isolation);

   no_host = rmfield(meta, 'hostname');
   returned = icemodel.test.helpers.perfMeasurementQuality(rows, no_host);
   testCase.verifyFalse(returned.conditions.one_machine);

   drifted = meta;
   drifted.ambient_stable = false;
   returned = icemodel.test.helpers.perfMeasurementQuality(rows, drifted);
   testCase.verifyFalse(returned.conditions.ambient_stable);
   testCase.verifyEqual(numel(returned.reasons), 1);

   overridden = meta;
   overridden.ambient_drift_accepted = true;
   returned = icemodel.test.helpers.perfMeasurementQuality(rows, overridden);
   testCase.verifyFalse(returned.conditions.ambient_stable);

   no_attestation = rmfield(meta, 'attestation');
   returned = icemodel.test.helpers.perfMeasurementQuality(rows, no_attestation);
   testCase.verifyFalse(returned.conditions.attestation);

   % An attestation without the run-start load, the shape a file written
   % before that field existed would carry, does not satisfy the condition.
   old_shape = meta;
   old_shape.attestation = rmfield(meta.attestation, 'load_average_at_start');
   returned = icemodel.test.helpers.perfMeasurementQuality(rows, old_shape);
   testCase.verifyFalse(returned.conditions.attestation);

   % A run that skipped the anchor names that cause, not a drift.
   unmeasured = meta;
   unmeasured.anchor_measured = false;
   unmeasured.ambient_stable = false;
   returned = icemodel.test.helpers.perfMeasurementQuality(rows, unmeasured);
   testCase.verifyFalse(returned.conditions.ambient_stable);
   testCase.verifyEqual(returned.reasons, "the ambient anchor was not measured");
end

function test_quality_ignores_the_comparison_verdict(testCase)
   % A run whose comparison failed still has trustworthy measurements, and
   % a drifted run keeps every case's own verdict and numbers.
   [rows, meta] = compliantSource();
   rows.passed_perf = [false; true];
   returned = icemodel.test.helpers.perfMeasurementQuality(rows, meta);
   testCase.verifyTrue(returned.passed);

   drifted = meta;
   drifted.ambient_stable = false;
   returned = icemodel.test.helpers.perfMeasurementQuality(rows, drifted);
   testCase.verifyFalse(returned.passed);
   testCase.verifyEqual(rows.passed_perf, [false; true]);
   testCase.verifyEqual(rows.median_wall_s, [10; 20]);
end

function test_release_source_check_accepts_a_compliant_source(testCase)
   % The compliant source passes on its own machine.
   [rows, meta] = compliantSource();
   testCase.verifyWarningFree(@() ...
      icemodel.test.helpers.assertReleasePerfBaselineSource(rows, meta, ...
      current_identity="macbook-air-2"));
   % A legacy .local hostname normalizes to the same identity.
   legacy = meta;
   legacy.hostname = "MacBook-Air-2.local";
   testCase.verifyWarningFree(@() ...
      icemodel.test.helpers.assertReleasePerfBaselineSource(rows, legacy, ...
      current_identity="macbook-air-2"));
end

function test_release_source_check_refuses_each_condition(testCase)
   % Every release-only condition refuses on its own, with its own error id.
   [rows, meta] = compliantSource();
   check = @(r, m) icemodel.test.helpers.assertReleasePerfBaselineSource( ...
      r, m, current_identity="macbook-air-2");

   overridden = meta;
   overridden.ambient_drift_accepted = true;
   testCase.verifyError(@() check(rows, overridden), ...
      'icemodel:test:releasePerfSourceQuality');

   % Self-consistent, but measured on another machine.
   foreign_host = meta;
   foreign_host.hostname = "other-machine";
   testCase.verifyError(@() check(rows, foreign_host), ...
      'icemodel:test:releasePerfSourceForeignHost');

   busy = meta;
   busy.attestation = icemodel.test.helpers.summarizeMachineState( ...
      fakeSamples([1.0 1.0], [0 1], [true true]));
   testCase.verifyError(@() check(rows, busy), ...
      'icemodel:test:releasePerfSourceForeignMatlab');

   % The run-start load is gated; the run's own subprocesses raise the
   % later samples, which are recorded but not gated.
   loaded = meta;
   loaded.attestation = icemodel.test.helpers.summarizeMachineState( ...
      fakeSamples([4.5 1.0], [0 0], [true true]));
   testCase.verifyError(@() check(rows, loaded), ...
      'icemodel:test:releasePerfSourceLoadAverage');
   busy_later = meta;
   busy_later.attestation = icemodel.test.helpers.summarizeMachineState( ...
      fakeSamples([1.5 6.4 5.8], [0 0 0], [true true true]));
   testCase.verifyWarningFree(@() check(rows, busy_later));

   battery = meta;
   battery.attestation = icemodel.test.helpers.summarizeMachineState( ...
      fakeSamples([1.0 1.0], [0 0], [true false]));
   testCase.verifyError(@() check(rows, battery), ...
      'icemodel:test:releasePerfSourceBatteryPower');

   % A NaN attestation value, from a failed probe, is not a quiet machine.
   unknown = meta;
   unknown.attestation = icemodel.test.helpers.summarizeMachineState( ...
      fakeSamples([1.0 NaN], [0 NaN], [true true]));
   testCase.verifyError(@() check(rows, unknown), ...
      'icemodel:test:releasePerfSourceForeignMatlab');
end

function test_release_load_threshold_lives_in_the_policy(testCase)
   % The load-average ceiling is read from perfMeasurementPolicy only.
   policy = icemodel.test.helpers.perfMeasurementPolicy();
   testCase.verifyEqual(policy.load_average_max, 4.0);
   [rows, meta] = compliantSource();
   meta.attestation = icemodel.test.helpers.summarizeMachineState( ...
      fakeSamples([4.0 1.0], [0 0], [true true]));
   testCase.verifyWarningFree(@() ...
      icemodel.test.helpers.assertReleasePerfBaselineSource(rows, meta, ...
      current_identity="macbook-air-2"));
   strict = policy;
   strict.load_average_max = 3.0;
   testCase.verifyError(@() ...
      icemodel.test.helpers.assertReleasePerfBaselineSource(rows, meta, ...
      current_identity="macbook-air-2", policy=strict), ...
      'icemodel:test:releasePerfSourceLoadAverage');
end

function test_managed_perf_snapshot_refuses_a_noncompliant_source(testCase)
   % snapshotBaseline applies the release-source check on the managed path
   % before it writes anything, so a rolling file with an accepted drift
   % override cannot become a release file.
   fixture = testCase.applyFixture( ...
      matlab.unittest.fixtures.TemporaryFolderFixture);
   source_file = fullfile(fixture.Folder, 'perf_rolling.mat');
   [~, meta] = compliantSource();
   % The forcing check runs before the release-source check, so the rows
   % must carry the rolling forcing identity of the registered release.
   PerfBaseline = icemodel.test.helpers.loadBaseline("perf", ...
      smbmodel="icemodel", simyear=2016);
   meta.git_revision = "v1.0.0-600-gabcdef12";
   meta.ambient_drift_accepted = true;
   save(source_file, 'PerfBaseline', 'meta');
   testCase.verifyError(@() icemodel.test.helpers.snapshotBaseline( ...
      "perf", "v1.3", "icemodel", false, string.empty(), 2016, ...
      source_file), 'icemodel:test:releasePerfSourceQuality');
   managed_file = icemodel.test.helpers.baselineFilePath( ...
      "perf", smbmodel="icemodel", baseline_type="release", ...
      baseline_tag="v1.3", simyear=2016);
   testCase.verifyFalse(isfile(managed_file));

   % Spelling the managed path out as output_file does not skip the gate.
   testCase.verifyError(@() icemodel.test.helpers.snapshotBaseline( ...
      "perf", "v1.3", "icemodel", false, string(managed_file), 2016, ...
      source_file), 'icemodel:test:releasePerfSourceQuality');
   testCase.verifyFalse(isfile(managed_file));
end

function test_build_quality_gate_allows_only_the_drift_override(testCase)
   % A managed build refuses any failed quality condition except a drifted
   % anchor whose drift the caller accepted.
   [rows, meta] = compliantSource();
   quality = icemodel.test.helpers.perfMeasurementQuality(rows, meta);
   testCase.verifyWarningFree(@() ...
      icemodel.test.helpers.assertPerfBuildQuality(quality, false));

   drifted = meta;
   drifted.ambient_stable = false;
   quality = icemodel.test.helpers.perfMeasurementQuality(rows, drifted);
   testCase.verifyError(@() ...
      icemodel.test.helpers.assertPerfBuildQuality(quality, false), ...
      'icemodel:test:perf:buildQuality');
   testCase.verifyWarningFree(@() ...
      icemodel.test.helpers.assertPerfBuildQuality(quality, true));

   session = drifted;
   session.isolation = "session";
   quality = icemodel.test.helpers.perfMeasurementQuality(rows, session);
   testCase.verifyError(@() ...
      icemodel.test.helpers.assertPerfBuildQuality(quality, true), ...
      'icemodel:test:perf:buildQuality');

   invalid_rows = rows;
   invalid_rows.valid(1) = false;
   quality = icemodel.test.helpers.perfMeasurementQuality(invalid_rows, meta);
   testCase.verifyError(@() ...
      icemodel.test.helpers.assertPerfBuildQuality(quality, false), ...
      'icemodel:test:perf:buildQuality');
end

function test_managed_session_build_refuses_before_measuring(testCase)
   % A session build of the managed rolling file refuses at the entry point,
   % before any bootstrap or model run; a custom output file may use it.
   testCase.verifyError(@() build_perf_baseline( ...
      isolation="session", include_benchmarks=false), ...
      'icemodel:test:perf:managedBuildRequiresProcessIsolation');

   % Naming the managed rolling path as output_file is still a managed build.
   managed_file = icemodel.test.helpers.baselineFilePath("perf", ...
      smbmodel="icemodel", simyear=2016);
   testCase.verifyError(@() build_perf_baseline( ...
      isolation="session", include_benchmarks=false, ...
      smbmodel="icemodel", output_file=string(managed_file)), ...
      'icemodel:test:perf:managedBuildRequiresProcessIsolation');
end
