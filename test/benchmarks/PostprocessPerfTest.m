classdef PostprocessPerfTest < matlab.perftest.TestCase
   %POSTPROCESSPERFTEST Benchmark postprocess hotspots on realistic output.

   properties
      workspace
      ice1_tt
   end

   methods (TestClassSetup)
      function buildPostprocessInput(testCase)
         % Run one short synthetic quarter-hour case, then tile its surface
         % output into a full-year-sized timetable for the retime benchmarks.
         testCase.workspace = icemodel.test.fixtures.makeSyntheticWorkspace( ...
            2016, configure=true, nsteps=96, dt_seconds=900);
         opts = icemodel.test.helpers.buildSyntheticOpts( ...
            testCase.workspace, 'skinmodel', 2016, dt=900, solver=1, ...
            testname='postprocess_perf_kernel');
         [ice1_raw, ~, opts] = icemodel.test.helpers.runSmbModel(opts);
         met = icemodel.loadmet(opts);
         ice1_source = rawIce1ToTimetable(ice1_raw, met.Time);
         testCase.ice1_tt = repeatQuarterHourTimetable(ice1_source, 35040);
      end
   end

   methods (TestClassTeardown)
      function cleanupPostprocessInput(testCase)
         % Tear down the synthetic workspace after the class finishes.
         icemodel.test.fixtures.cleanupSyntheticWorkspace(testCase.workspace);
      end
   end

   methods (Test)
      function testHourlyRetimeLegacy(testCase)
         % Benchmark the timetable RETIME branch used historically.
         batch_size = 4;
         ice1_hourly = legacyHourlyMean(testCase.ice1_tt);
         testCase.assertTrue(~isempty(ice1_hourly));

         while testCase.keepMeasuring
            for n = 1:batch_size
               ice1_hourly = legacyHourlyMean(testCase.ice1_tt);
            end
            if isempty(ice1_hourly)
               error('legacy RETIME benchmark produced an empty timetable')
            end
         end
      end

      function testHourlyRetimeFixedStep(testCase)
         % Benchmark the fixed-step hourly aggregation helper.
         batch_size = 64;
         ice1_hourly = icemodel.internal.retimeHourlyFixedStep(testCase.ice1_tt);
         testCase.assertTrue(~isempty(ice1_hourly));

         while testCase.keepMeasuring
            for n = 1:batch_size
               ice1_hourly = icemodel.internal.retimeHourlyFixedStep( ...
                  testCase.ice1_tt);
            end
            if isempty(ice1_hourly)
               error(['fixed-step hourly benchmark produced an empty ', ...
                  'timetable'])
            end
         end
      end
   end
end

function ice1_tt = rawIce1ToTimetable(ice1_raw, time)
   %RAWICE1TOTIMETABLE Convert raw model output to POSTPROCESS timetable form.

   % Match POSTPROCESS by converting logical convergence flags to single
   % before the timetable conversion.
   if isfield(ice1_raw, 'Tice_converged')
      ice1_raw.Tice_converged = single(ice1_raw.Tice_converged);
   end
   if isfield(ice1_raw, 'Tsfc_converged')
      ice1_raw.Tsfc_converged = single(ice1_raw.Tsfc_converged);
   end

   % Build the timetable shape timed by the hourly aggregation hotspot.
   time.TimeZone = 'UTC';
   ice1_tt = struct2table(ice1_raw);
   ice1_tt = table2timetable(ice1_tt, 'RowTimes', time);
end

function TT = repeatQuarterHourTimetable(TT, n_rows)
   %REPEATQUARTERHOURTIMETABLE Tile a short timetable to a target row count.

   % Repeat the source data enough times to hit the requested representative
   % benchmark length, then assign one continuous quarter-hour time axis.
   vars = TT.Properties.VariableNames;
   data = TT{:, vars};
   n_repeat = ceil(n_rows / size(data, 1));
   data = repmat(data, n_repeat, 1);
   data = data(1:n_rows, :);

   time = datetime(2016, 1, 1, 0, 0, 0, 'TimeZone', 'UTC') + ...
      minutes(15) * transpose(0:n_rows-1);
   TT = array2timetable(data, 'RowTimes', time, 'VariableNames', vars);
end

function ice1_hourly = legacyHourlyMean(ice1_tt)
   %LEGACYHOURLYMEAN Reproduce the timetable RETIME branch from POSTPROCESS.

   % Match the legacy path exactly, including the leap-day removal.
   ice1_hourly = retime(ice1_tt, 'hourly', 'mean');
   ice1_hourly = ice1_hourly(~(month(ice1_hourly.Time) == 2 ...
      & day(ice1_hourly.Time) == 29), :);
end
