function tests = test_forcing_batch_equivalence
   %TEST_FORCING_BATCH_EQUIVALENCE Multi-point builds equal single-point loops.
   %
   % The 1ps.19 performance refactor lets the gridded forcing builders
   % (buildMarData/buildMerraData/buildRacmoData and the buildMar/MerraMet
   % wrappers) accept a LIST of points (Nx2 [lat lon]) and open each source
   % file ONCE for the whole list, rather than re-opening it per point. This
   % suite is the hard equivalence gate: building 2-3 points the NEW way (one
   % multi-point call) must be byte-identical (isequaln on the variables,
   % times, and CustomProperties) to the OLD way (a loop of single-point
   % calls). Payload metadata is part of that equality contract. Lanes self-skip
   % when their staged fast fixtures are not on disk.
   %
   % Reads the small fixture subset under test/data/forcing. This belongs in
   % regression because it exercises real RCM I/O and batch-vs-loop behavior,
   % not an isolated unit.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % A small cluster of points near the KAN transect (western Greenland
   % ablation/percolation), so every point resolves on every grid.
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;

   testCase.TestData.points = [ ...
      67.0955, -50.0699;    % ~KAN_L
      67.0670, -48.8355;    % ~KAN_M
      67.0003, -47.0245];   % ~KAN_U
   testCase.TestData.year = 2012;
   cfg = icemodel.config('getenv', true);
   forcing_root = string(fullfile(cfg.ICEMODEL_DATA_PATH, 'forcing'));

   testCase.TestData.mar = firstWithData( ...
      fullfile(forcing_root, 'mar'), ...
      @(p) ~isempty(dir(fullfile(p, "MARv3.11*.nc"))));
   testCase.TestData.merra = firstWithData( ...
      fullfile(forcing_root, 'merra2'), ...
      @(p) ~isempty(dir(fullfile(p, "slv", "*_Nx.*.nc4*"))));
   testCase.TestData.racmo = firstWithData( ...
      fullfile(forcing_root, 'racmo'), ...
      @(p) ~isempty(dir(fullfile(p, "*.RACMO23p3_*.nc"))));
end

function test_mar_batch_equals_single_loop(testCase)
   % buildMarData/buildMarMet: a 3-point list equals the single-point loop.
   src = testCase.TestData.mar;
   testCase.assumeTrue(strlength(src) > 0, ...
      'MAR fixture data not available under test/data/forcing');
   pts = testCase.TestData.points;
   yr = testCase.TestData.year;

   batch = icemodel.forcing.buildMarData(pts, yr, source_dir=src);
   verifyBatchEqualsLoop(testCase, batch, pts, ...
      @(p) icemodel.forcing.buildMarData(p, yr, source_dir=src), 'MAR Data');

   % The thin met wrapper threads the point list through data2met too.
   [metbatch, metmetadata, metData] = ...
      icemodel.forcing.buildMarMet(pts, yr, source_dir=src);
   verifyBatchEqualsLoop(testCase, metbatch, pts, ...
      @(p) icemodel.forcing.buildMarMet(p, yr, source_dir=src), 'MAR met');
   verifyMetBatchOutputs(testCase, metbatch, metmetadata, metData, batch, ...
      'MAR met');
end

function test_merra_batch_equals_single_loop(testCase)
   % buildMerraData/buildMerraMet: a 3-point list equals the single-point loop.
   src = testCase.TestData.merra;
   testCase.assumeTrue(strlength(src) > 0, ...
      'MERRA-2 fixture data not available under test/data/forcing');
   pts = testCase.TestData.points;
   yr = testCase.TestData.year;

   batch = icemodel.forcing.buildMerraData(pts, yr, source_dir=src);
   verifyBatchEqualsLoop(testCase, batch, pts, ...
      @(p) icemodel.forcing.buildMerraData(p, yr, source_dir=src), 'MERRA Data');

   [metbatch, metmetadata, metData] = ...
      icemodel.forcing.buildMerraMet(pts, yr, source_dir=src);
   verifyBatchEqualsLoop(testCase, metbatch, pts, ...
      @(p) icemodel.forcing.buildMerraMet(p, yr, source_dir=src), 'MERRA met');
   verifyMetBatchOutputs(testCase, metbatch, metmetadata, metData, batch, ...
      'MERRA met');
end

function test_racmo_batch_equals_single_loop(testCase)
   % buildRacmoData: a 3-point list equals the single-point loop.
   src = testCase.TestData.racmo;
   testCase.assumeTrue(strlength(src) > 0, ...
      'RACMO fixture data not available under test/data/forcing');
   pts = testCase.TestData.points;
   yr = testCase.TestData.year;

   batch = icemodel.forcing.buildRacmoData(pts, yr, source_dir=src);
   verifyBatchEqualsLoop(testCase, batch, pts, ...
      @(p) icemodel.forcing.buildRacmoData(p, yr, source_dir=src), 'RACMO Data');
end

function test_single_point_still_returns_a_timetable(testCase)
   % N=1 is the single-point path: a 1x2 point returns ONE timetable (not a
   % 1x1 cell), preserving the legacy single-location contract.
   src = testCase.TestData.mar;
   testCase.assumeTrue(strlength(src) > 0, ...
      'MAR fixture data not available under test/data/forcing');
   one = icemodel.forcing.buildMarData( ...
      testCase.TestData.points(1, :), testCase.TestData.year, source_dir=src);
   testCase.verifyTrue(istimetable(one));
end

%% Local helpers
function verifyBatchEqualsLoop(testCase, batch, pts, buildOne, label)
   %VERIFYBATCHEQUALSLOOP Assert a 1xN cell batch equals N single-point builds.
   testCase.verifyTrue(iscell(batch), ...
      sprintf('%s: multi-point build must return a cell array', label));
   testCase.verifyEqual(numel(batch), size(pts, 1), ...
      sprintf('%s: one Data per point', label));
   for p = 1:size(pts, 1)
      single = buildOne(pts(p, :));
      verifyTimetablesIdentical(testCase, batch{p}, single, ...
         sprintf('%s point %d', label, p));
   end
end

function verifyMetBatchOutputs(testCase, met, metadata, Data, expectedData, label)
   %VERIFYMETBATCHOUTPUTS Enforce the common three-output met-builder contract.
   testCase.verifyTrue(isstruct(metadata), ...
      sprintf('%s: metadata must be a struct array', label));
   testCase.verifySize(metadata, size(met), ...
      sprintf('%s: metadata must preserve the met collection shape', label));
   testCase.verifySize(Data, size(met), ...
      sprintf('%s: Data must preserve the met collection shape', label));
   for k = 1:numel(met)
      testCase.verifyEqual(metadata(k), met{k}.Properties.UserData, ...
         sprintf('%s: metadata %d must equal met UserData', label, k));
      verifyTimetablesIdentical(testCase, Data{k}, expectedData{k}, ...
         sprintf('%s source Data %d', label, k));
   end
end

function verifyTimetablesIdentical(testCase, a, b, label)
   %VERIFYTIMETABLESIDENTICAL isequaln on variables, times, units, and
   % the userdata metadata/CustomProperties.
   testCase.verifyEqual(string(a.Properties.VariableNames), ...
      string(b.Properties.VariableNames), ...
      sprintf('%s: variable names differ', label));
   testCase.verifyTrue(isequaln(a.Time, b.Time), ...
      sprintf('%s: time axes differ', label));
   testCase.verifyTrue(isequaln(timetable2table(a), timetable2table(b)), ...
      sprintf('%s: variable data differ', label));
   testCase.verifyEqual(string(a.Properties.VariableUnits), ...
      string(b.Properties.VariableUnits), ...
      sprintf('%s: variable units differ', label));
   testCase.verifyTrue(isequaln(a.Properties.UserData, ...
      b.Properties.UserData), ...
      sprintf('%s: UserData metadata differ', label));

   props = ["X", "Y", "Lat", "Lon", "Elev", "Slope", "ScalarUnits"];
   for q = props
      if isprop(a.Properties.CustomProperties, q)
         testCase.verifyTrue(isequaln( ...
            a.Properties.CustomProperties.(q), ...
            b.Properties.CustomProperties.(q)), ...
            sprintf('%s: CustomProperties.%s differ', label, q));
      end
   end
end

function p = firstWithData(candidates, hasData)
   %FIRSTWITHDATA First candidate dir that exists and passes the data probe.
   p = "";
   for c = candidates
      if isfolder(c) && hasData(c)
         p = c;
         return
      end
   end
end
