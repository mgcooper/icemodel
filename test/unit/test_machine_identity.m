function tests = test_machine_identity
   %TEST_MACHINE_IDENTITY Verify the network-independent machine identity.
   %
   % icemodel.test.helpers.machineHostname no longer trusts the
   % network-dependent `hostname` command on macOS: it prefers
   % `scutil --get LocalHostName`, and
   % icemodel.test.helpers.normalizeMachineIdentity folds any raw name
   % (from either probe, or a saved baseline's meta.hostname) to one
   % stable, comparable identity (DesignSpec decisions 18 and 19). These
   % tests cover every machineHostname probe branch, the normalization
   % rules, perfBaselineCompatibility's use of the normalized identity on
   % both sides of a comparison, and run_aa_acceptance's normalized
   % hostname uniqueness check against synthetic artifacts.
   %
   % See also: icemodel.test.helpers.machineHostname,
   %  icemodel.test.helpers.normalizeMachineIdentity,
   %  icemodel.test.helpers.perfBaselineCompatibility,
   %  icemodel.test.helpers.formalPerformanceVerdict, run_aa_acceptance
   tests = functiontests(localfunctions);
end

%% normalizeMachineIdentity

function test_normalize_strips_legacy_local_suffix(testCase)
   % The legacy `hostname`-derived value must normalize to the same
   % identity as the new scutil-derived value.
   returned = icemodel.test.helpers.normalizeMachineIdentity( ...
      "MacBook-Air-2.local");
   expected = "macbook-air-2";
   testCase.verifyEqual(returned, expected);
end

function test_normalize_leaves_a_bare_name_unchanged(testCase)
   % A name with no case, whitespace, or suffix to fold is returned as-is.
   returned = icemodel.test.helpers.normalizeMachineIdentity("x");
   expected = "x";
   testCase.verifyEqual(returned, expected);
end

function test_normalize_is_case_insensitive(testCase)
   % An uppercase raw name and its ".LOCAL" suffix both fold away.
   returned = icemodel.test.helpers.normalizeMachineIdentity("X.LOCAL");
   expected = "x";
   testCase.verifyEqual(returned, expected);
end

function test_normalize_strips_one_trailing_suffix_only(testCase)
   % Only the trailing ".local" is a Bonjour/mDNS suffix; an interior one
   % is part of the name and stays.
   returned = icemodel.test.helpers.normalizeMachineIdentity( ...
      "a.local.local");
   expected = "a.local";
   testCase.verifyEqual(returned, expected);
end

function test_normalize_trims_surrounding_whitespace(testCase)
   % A probe result padded with whitespace (as a shell command returns)
   % still normalizes cleanly.
   returned = icemodel.test.helpers.normalizeMachineIdentity( ...
      sprintf('  MacBook-Air-2.local  \n'));
   expected = "macbook-air-2";
   testCase.verifyEqual(returned, expected);
end

function test_normalize_blank_input_stays_blank(testCase)
   % A blank value carries no identity; callers still see it as blank
   % through their own isblanktext check.
   returned = icemodel.test.helpers.normalizeMachineIdentity("   ");
   testCase.verifyTrue(isblanktext(returned));
end

%% machineHostname probe branches

function test_hostname_prefers_scutil_on_macos(testCase)
   % macOS reads scutil's LocalHostName first. The fake probe used here
   % errors if machineHostname falls back to `hostname`, so a passing
   % test proves the fallback is not taken when scutil succeeds.
   returned = icemodel.test.helpers.machineHostname( ...
      @fakeProbeScutilOnly, true);
   expected = "macbook-air-2";
   testCase.verifyEqual(returned, expected);
end

function test_hostname_falls_back_when_scutil_fails(testCase)
   % A failed scutil probe on macOS falls back to `hostname`; the legacy
   % ".local"-suffixed value it returns still normalizes to the same
   % identity scutil would have produced.
   returned = icemodel.test.helpers.machineHostname( ...
      @fakeProbeScutilFailsHostnameSucceeds, true);
   expected = "macbook-air-2";
   testCase.verifyEqual(returned, expected);
end

function test_hostname_uses_only_hostname_off_macos(testCase)
   % A non-macOS platform never probes scutil: is_mac=false skips that
   % branch in machineHostname entirely, so the fake probe's scutil-error
   % path is never reached.
   returned = icemodel.test.helpers.machineHostname( ...
      @fakeProbeHostnameOnlyErrorsOnScutil, false);
   expected = "some-host";
   testCase.verifyEqual(returned, expected);
end

function test_hostname_rejects_command_failure_on_both_probes(testCase)
   % Neither probe can identify the machine when both fail.
   testCase.verifyError(@() icemodel.test.helpers.machineHostname( ...
      @fakeProbeBothFail, true), ...
      'icemodel:test:perf:machineIdentityUnavailable');
end

function test_hostname_rejects_blank_output_on_both_probes(testCase)
   % A successful command with no machine name cannot identify the host,
   % whether it comes from scutil or the hostname fallback.
   testCase.verifyError(@() icemodel.test.helpers.machineHostname( ...
      @fakeProbeBothBlank, true), ...
      'icemodel:test:perf:machineIdentityUnavailable');
end

%% perfBaselineCompatibility with the normalized identity

function test_compatibility_survives_a_simulated_network_change(testCase)
   % A baseline saved under the legacy `hostname` value stays compatible
   % with a current identity produced by the new LocalHostName-derived
   % probe: the simulated-network-change case and the legacy-value case
   % that DesignSpec decisions 18-19 require to keep pre-change rolling
   % baselines comparable.
   meta = struct('matlab_version', string(version), ...
      'host', string(computer), 'hostname', "MacBook-Air-2.local", ...
      'isolation', "process");
   [returned, reason] = icemodel.test.helpers.perfBaselineCompatibility( ...
      meta, "process", "macbook-air-2");
   testCase.verifyTrue(returned);
   testCase.verifyEqual(reason, "");
end

function test_compatibility_normalizes_an_injected_current_identity(testCase)
   % A caller-supplied CURRENT_IDENTITY is normalized the same way as the
   % live probe, so a raw, un-normalized value still compares correctly.
   meta = struct('matlab_version', string(version), ...
      'host', string(computer), 'hostname', "macbook-air-2", ...
      'isolation', "process");
   [returned, reason] = icemodel.test.helpers.perfBaselineCompatibility( ...
      meta, "process", "  MacBook-Air-2.LOCAL  ");
   testCase.verifyTrue(returned);
   testCase.verifyEqual(reason, "");
end

function test_compatibility_rejects_a_genuinely_different_machine(testCase)
   % Two distinct machines are not comparable, and the reason names both,
   % using their normalized identities.
   meta = struct('matlab_version', string(version), ...
      'host', string(computer), 'hostname', "workstation-one", ...
      'isolation', "process");
   [returned, reason] = icemodel.test.helpers.perfBaselineCompatibility( ...
      meta, "process", "workstation-two");
   testCase.verifyFalse(returned);
   testCase.verifySubstring(reason, "workstation-one");
   testCase.verifySubstring(reason, "workstation-two");
end

%% DesignSpec decision 19: an incompatible host still writes the record

function test_incompatible_host_annotates_but_does_not_block_recording(testCase)
   % perfBaselineCompatibility marks the row incompatible and states why.
   % test/run_perf_suite.m stores meta.baseline_compatible and
   % meta.compare_reason unconditionally after this call (see its meta
   % assignments around line 520), and it saves the measured
   % median_wall_s for every case regardless of the verdict: an
   % incompatible host annotates the record, it does not skip writing it.
   % formalPerformanceVerdict mirrors that: it fails closed but never
   % erases the caller's own measured value.
   meta = struct('matlab_version', string(version), ...
      'host', string(computer), 'hostname', "workstation-one", ...
      'isolation', "process");
   [baseline_compatible, compare_reason] = ...
      icemodel.test.helpers.perfBaselineCompatibility( ...
      meta, "process", "workstation-two");
   testCase.verifyFalse(baseline_compatible);

   current_wall = 12.5;
   [passed, ref_wall, floor_wall, gate_wall, reason] = ...
      icemodel.test.helpers.formalPerformanceVerdict( ...
      true, current_wall, table(), [], baseline_compatible, 0.2, ...
      compare_reason);

   % The verdict fails closed and carries no reference or gate (there is
   % no compatible baseline row to compare against), but the caller's own
   % measured CURRENT_WALL is untouched and still available to record.
   % formalPerformanceVerdict converts its REASON to string(); compare
   % against the same type rather than perfBaselineCompatibility's char.
   testCase.verifyFalse(passed);
   testCase.verifyEqual(reason, string(compare_reason));
   testCase.verifyTrue(isnan(ref_wall) && isnan(floor_wall) ...
      && isnan(gate_wall));
   testCase.verifyEqual(current_wall, 12.5);
end

%% run_aa_acceptance: normalized hostname uniqueness on synthetic artifacts

function test_aa_acceptance_accepts_legacy_and_normalized_hostnames(testCase)
   % Two synthetic artifacts recorded with the legacy hostname value and
   % the new LocalHostName-derived identity must pass the hostname
   % uniqueness check: the "legacy value" case DesignSpec decisions 18-19
   % require, so pre-change rolling artifacts stay comparable with
   % artifacts saved after the probe change. run_aa_acceptance loads
   % artifacts from disk, so the pair is written under this test's own
   % temporary folder fixture rather than compared in memory.
   fixture = testCase.applyFixture( ...
      matlab.unittest.fixtures.TemporaryFolderFixture);

   meta = identityArtifactMeta("run-a", "MacBook-Air-2.local");
   case_summary = identityCaseSummary(100.0, true);
   file_a = fullfile(fixture.Folder, 'identity_a.mat');
   save(file_a, 'meta', 'case_summary');

   meta = identityArtifactMeta("run-b", "macbook-air-2");
   case_summary.median_wall_s = 100.2;
   file_b = fullfile(fixture.Folder, 'identity_b.mat');
   save(file_b, 'meta', 'case_summary');

   % Discard the printed verdict so the suite log stays quiet.
   [~] = evalc('report = run_aa_acceptance(file_a, file_b);');
   testCase.verifyTrue(report.passed);
end

function test_aa_acceptance_rejects_a_genuinely_different_machine(testCase)
   % Two artifacts from distinct machines must still fail the hostname
   % uniqueness check after normalization; normalization must not paper
   % over a real environment mismatch.
   fixture = testCase.applyFixture( ...
      matlab.unittest.fixtures.TemporaryFolderFixture);

   meta = identityArtifactMeta("run-a", "workstation-one");
   case_summary = identityCaseSummary(100.0, true);
   file_a = fullfile(fixture.Folder, 'identity_a.mat');
   save(file_a, 'meta', 'case_summary');

   meta = identityArtifactMeta("run-b", "workstation-two");
   case_summary.median_wall_s = 100.2;
   file_b = fullfile(fixture.Folder, 'identity_b.mat');
   save(file_b, 'meta', 'case_summary');

   testCase.verifyError(@() run_aa_acceptance(file_a, file_b), ...
      'icemodel:test:aaAcceptance:environmentMismatch');
end

%% Local fakes and fixtures

function [status, value] = fakeProbeScutilOnly(command)
   %FAKEPROBESCUTILONLY Answer scutil only; error if hostname is probed too.
   if contains(command, "scutil")
      status = 0;
      value = "MacBook-Air-2";
   else
      error('icemodel:test:machineIdentity:unexpectedProbe', ...
         ['machineHostname must not fall back to hostname when ' ...
         'scutil succeeds']);
   end
end

function [status, value] = fakeProbeScutilFailsHostnameSucceeds(command)
   %FAKEPROBESCUTILFAILSHOSTNAMESUCCEEDS Fail scutil, succeed on hostname.
   if contains(command, "scutil")
      status = 1;
      value = "";
   else
      status = 0;
      value = "MacBook-Air-2.local";
   end
end

function [status, value] = fakeProbeHostnameOnlyErrorsOnScutil(command)
   %FAKEPROBEHOSTNAMEONLYERRORSONSCUTIL Prove the non-macOS path skips scutil.
   if contains(command, "scutil")
      error('icemodel:test:machineIdentity:unexpectedProbe', ...
         'machineHostname must not probe scutil off macOS');
   end
   status = 0;
   value = "some-host";
end

function [status, value] = fakeProbeBothFail(command)
   %FAKEPROBEBOTHFAIL Fail every probe call, whatever COMMAND asks for.
   assert(isstring(command) || ischar(command));
   status = 1;
   value = "";
end

function [status, value] = fakeProbeBothBlank(command)
   %FAKEPROBEBOTHBLANK Succeed but report no name, whatever COMMAND asks for.
   assert(isstring(command) || ischar(command));
   status = 0;
   value = "   ";
end

function meta = identityArtifactMeta(run_name, hostname)
   %IDENTITYARTIFACTMETA Build the metadata run_aa_acceptance requires.
   %
   % HOSTNAME is left exactly as supplied (not pre-normalized), so the
   % caller's test exercises run_aa_acceptance's own normalization.
   meta = struct('ambient_stable', true, 'isolation', "process", ...
      'run_name', string(run_name), 'hostname', string(hostname), ...
      'matlab_version', string(version), 'git_revision', "test-rev", ...
      'tier', "smoke", 'simyear', 2016, 'n_runs', 3, 'n_warmups', 1, ...
      'tol_perf', 0.2, 'data_root', "test-data");
end

function case_summary = identityCaseSummary(median_wall_s, valid)
   %IDENTITYCASESUMMARY Build one complete A/A workload row.
   case_summary = table("icemodel_kanm_2016_solver1", ...
      "promice_filled", median_wall_s, valid, ...
      'VariableNames', {'case_id', 'forcings', 'median_wall_s', 'valid'});
end
