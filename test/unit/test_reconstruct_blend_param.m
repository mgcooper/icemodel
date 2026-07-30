classdef test_reconstruct_blend_param < matlab.unittest.TestCase
   %TEST_RECONSTRUCT_BLEND_PARAM Parameterized seam-blend coverage.
   %
   % Sweeps the seam-blend regimes the scalar-cased suite cannot cover
   % combinatorially: run lengths spanning the single-sample,
   % shorter-than-taper, boundary (exactly two taper windows), and
   % disjoint-ramp regimes, crossed with offset directions. Every
   % combination must fill with the offset method and land both anchored
   % seams inside the per-run jump limit.

   properties (TestParameter)
      % Run lengths in hours, chosen to hit each blend regime for the
      % default 6 h taper: 1 (single-sample), 3 and 8 (shorter than two
      % windows), 12 (exactly two windows), 48 (disjoint ramps).
      run_hours = struct('one', 1, 'three', 3, 'eight', 8, ...
         'twelve', 12, 'fortyeight', 48);
      % Method offset applied to the truth, in kelvin, both directions.
      offset = struct('plus', 25, 'minus', -25);
   end

   methods (TestClassSetup)
      function installPaths(testCase)
         % Install the verification path for namespace resolution.
         [~, ~, ~, ~, cleanup] = ...
            icemodel.test.helpers.bootstrapTestEnvironment();
         testCase.addTeardown(@() delete_handle(cleanup));
      end
   end

   methods (Test)
      function blended_fill_lands_inside_seam_limits(testCase, ...
            run_hours, offset)
         % An in-bounds estimate offset from the truth on both seams
         % must fill the whole run, keep the method's provenance code,
         % and land both anchored seams inside the per-run limit (the
         % base seasonal limit or the bridge floor, whichever is
         % larger). Tier 1 is disabled so every run length reaches the
         % method walk.
         series = icemodel.test.fixtures.makeReconstructSeries();
         truth = series.tair;
         codes = icemodel.forcing.reconstruct.provenanceCodes();
         n_total = numel(truth);
         gap = (4000:4000 + run_hours - 1).';
         series.tair(gap) = NaN;

         est = nan(n_total, 1);
         est(gap) = truth(gap) + offset;
         method = struct('name', "offset", ...
            'code', codes.station_transfer, 'estimate', est, ...
            'seasons', "all", 'buckets', []);
         cm = struct('channel', "tair", 'methods', method);

         result = icemodel.forcing.reconstruct.reconstructSeries( ...
            series, cm, interp_channels=string.empty(1, 0));

         returned = result.series.tair;
         testCase.verifyEqual(result.provenance.tair(gap), ...
            repmat(codes.station_transfer, numel(gap), 1));
         testCase.verifyTrue(all(isfinite(returned(gap))));

         % The seam limit the engine used: the seasonal base limit
         % floored at the linear bridge's per-step change.
         scale = icemodel.forcing.reconstruct.stepScale( ...
            series.Properties.RowTimes, series.tair);
         opts = icemodel.forcing.reconstruct.setopts();
         bridge = abs(truth(gap(end) + 1) - truth(gap(1) - 1)) ...
            / (numel(gap) + 1);
         limit = max(opts.jump_factor * scale.JJA, bridge * (1 + 1e-6));
         testCase.verifyLessThanOrEqual( ...
            abs(returned(gap(1)) - truth(gap(1) - 1)), limit);
         testCase.verifyLessThanOrEqual( ...
            abs(returned(gap(end)) - truth(gap(end) + 1)), limit);
      end
   end
end

function delete_handle(cleanup)
   %DELETE_HANDLE Destroy the environment cleanup object deterministically.
   delete(cleanup);
end
