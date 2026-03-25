%EXPLORE_RESULTS Demonstrate file path building and result loading helpers.
%
% This script shows how to use the test helpers to locate, load, and display
% baselines, artifacts, and references without running any test suite.
%
% Prerequisites:
%   Run icemodel.test.helpers.bootstrapTestEnvironment() first, or use fully
%   qualified names (icemodel.test.helpers.*) as shown below.

import icemodel.test.helpers.*

%% Baseline file paths

% Default rolling perf baseline for icemodel
baselineFilePath("perf")

% Rolling perf baseline for skinmodel
baselineFilePath("perf", smbmodel="skinmodel")

% Rolling regression baseline
baselineFilePath("regression")

% Release baseline by tag (tag implies baseline_type="release")
baselineFilePath("perf", baseline_tag="v1.1")
baselineFilePath("regression", baseline_tag="v1.1", smbmodel="skinmodel")

% Latest release baseline (auto-resolved by scanning baselines directory)
baselineFilePath("perf", baseline_type="release")
baselineFilePath("regression", baseline_type="release")

%% Load baselines

% Load the rolling perf baseline for all models
[perf_bl, perf_meta] = loadBaseline("perf");

% Load the rolling regression baseline for all models
reg_bl = loadBaseline("regression");

% Load a release baseline
[release_bl, release_meta] = loadBaseline("perf", baseline_tag="v1.1");

% Load a specific model
skin_bl = loadBaseline("perf", smbmodel="skinmodel");

% Load from an explicit file path
f = baselineFilePath("perf", smbmodel="skinmodel");
[data, meta] = loadBaseline("perf", filename=f);

%% Display baseline summaries

% Display perf baseline summary (works with tables or file paths)
displayPerfSummary(perf_bl)

% Display from a file path
displayPerfSummary(baselineFilePath("perf"))

% Display regression baseline summary
displayRegressionSummary(reg_bl)

%% Artifact file paths

% Most recent perf artifact (scans test/artifacts/ for latest run)
artifactFilePath("perf")

% Most recent regression artifact
artifactFilePath("regression")

% Artifact for a specific model and solver
artifactFilePath("perf", smbmodel="skinmodel", solver=1)

%% Load artifacts

% Load the most recent perf artifact
[perf_art, perf_art_meta] = loadArtifact("perf");

% Load the most recent regression artifact
[reg_art, reg_art_meta] = loadArtifact("regression");

% Display artifact results
displayPerfSummary(perf_art.case_summary, perf_art.benchmark)
displayRegressionSummary(reg_art.report)

%% Reference file paths and loading

% Runoff reference path
referenceFilePath("runoff")

% Load the runoff reference
ref = loadReference("runoff");
