classdef RenameRoundTest < matlab.perftest.TestCase
   %RENAMEROUNDTEST Formal rename/round microbenchmarks.
   %
   % This benchmark consolidates the legacy rename/round timing work into one
   % formal perftest target that uses synthetic but schema-representative data.
   %
   % Legacy note:
   %  - older exploratory variants included a script-based benchmark and a
   %    separate `time_roundData` helper that used real `ice1/ice2` outputs.
   %  - this reconciled benchmark does not use real data; the real-data note
   %    is preserved in the helper comments below where it affected the code
   %    choice, especially for postprocess-style rounding.
   %  - rounding is still active in the current production path:
   %    `icemodel.postprocess` calls `roundData(...)`.

   properties
      ice1_rename
      oldvars
      newvars
      ice2_round
      ice2lookup
      rename_batch_size
      round_batch_size
   end

   methods (TestMethodSetup)
      function generateTestData(testCase)
         % Build representative rename and rounding payloads once per method.
         [testCase.ice1_rename, testCase.oldvars, testCase.newvars] = ...
            createRenameData();
         [testCase.ice2_round, testCase.ice2lookup] = createRoundData();

         % Rename is fast enough to need batching for stable timing.
         testCase.rename_batch_size = 250;

         % The round variants are larger, but a small batch still reduces
         % sampling noise without changing the relative comparison much.
         testCase.round_batch_size = 2;
      end
   end

   methods (Test)
      function testRenameIsmember(testCase)
         % Measure the `ismember`-based rename strategy.
         renamed = renameUsingIsmember(testCase.ice1_rename, ...
            testCase.oldvars, testCase.newvars);
         testCase.assertNotEmpty(renamed)

         while testCase.keepMeasuring
            for n = 1:testCase.rename_batch_size
               renamed = renameUsingIsmember(testCase.ice1_rename, ...
                  testCase.oldvars, testCase.newvars);
            end
         end

         testCase.assertNotEmpty(renamed)
      end

      function testRenameIntersect(testCase)
         % Measure the `intersect`-based rename strategy.
         renamed = renameUsingIntersect(testCase.ice1_rename, ...
            testCase.oldvars, testCase.newvars);
         testCase.assertNotEmpty(renamed)

         while testCase.keepMeasuring
            for n = 1:testCase.rename_batch_size
               renamed = renameUsingIntersect(testCase.ice1_rename, ...
                  testCase.oldvars, testCase.newvars);
            end
         end

         testCase.assertNotEmpty(renamed)
      end

      function testRoundSwitch(testCase)
         % Measure the switch-based field-precision rounder.
         rounded = roundUsingSwitch(testCase.ice2_round);
         testCase.assertNotEmpty(rounded)

         while testCase.keepMeasuring
            for n = 1:testCase.round_batch_size
               rounded = roundUsingSwitch(testCase.ice2_round);
            end
         end

         testCase.assertNotEmpty(rounded)
      end

      function testRoundIsmember(testCase)
         % Measure the lookup-based rounder driven by `ismember`.
         rounded = roundUsingIsmember(testCase.ice2_round, ...
            testCase.ice2lookup);
         testCase.assertNotEmpty(rounded)

         while testCase.keepMeasuring
            for n = 1:testCase.round_batch_size
               rounded = roundUsingIsmember(testCase.ice2_round, ...
                  testCase.ice2lookup);
            end
         end

         testCase.assertNotEmpty(rounded)
      end

      function testRoundIntersect(testCase)
         % Measure the lookup-based rounder driven by `intersect`.
         rounded = roundUsingIntersect(testCase.ice2_round, ...
            testCase.ice2lookup);
         testCase.assertNotEmpty(rounded)

         while testCase.keepMeasuring
            for n = 1:testCase.round_batch_size
               rounded = roundUsingIntersect(testCase.ice2_round, ...
                  testCase.ice2lookup);
            end
         end

         testCase.assertNotEmpty(rounded)
      end

      function testRoundIsmemberPersistent(testCase)
         % Measure the persistent `ismember` lookup variant.
         rounded = roundUsingIsmemberPersistent(testCase.ice2_round);
         testCase.assertNotEmpty(rounded)

         while testCase.keepMeasuring
            for n = 1:testCase.round_batch_size
               rounded = roundUsingIsmemberPersistent(testCase.ice2_round);
            end
         end

         testCase.assertNotEmpty(rounded)
      end

      function testRoundIntersectPersistent(testCase)
         % Measure the persistent `intersect` lookup variant.
         rounded = roundUsingIntersectPersistent(testCase.ice2_round);
         testCase.assertNotEmpty(rounded)

         while testCase.keepMeasuring
            for n = 1:testCase.round_batch_size
               rounded = roundUsingIntersectPersistent(testCase.ice2_round);
            end
         end

         testCase.assertNotEmpty(rounded)
      end
   end
end

function [ice1, oldvars, newvars] = createRenameData()
   % Create a sample timetable for the rename benchmark.
   ice1 = timetable(datetime('now') + minutes(1:8760)', ...
      rand(8760, 1), rand(8760, 1), rand(8760, 1), rand(8760, 1), ...
      rand(8760, 1));
   ice1.Properties.VariableNames = {'Qsi', 'Qsr', 'Qsn', 'Qli', 'Qle'};

   oldvars = {'Qsi', 'Qsr', 'Qsn', 'Qli', 'Qle', 'Qln', 'Qh', 'Qe', ...
      'Qc', 'Qn', 'Tsfc'};
   newvars = {'swd', 'swu', 'swn', 'lwd', 'lwu', 'lwn', 'shf', 'lhf', ...
      'chf', 'netr', 'tsfc'};
end

function [ice2, ice2lookup] = createRoundData()
   % Create a sample structure for the rounding benchmark.
   ice2 = struct('f_ice', rand(500, 8760), 'Tice', rand(500, 8760), ...
      'cp_sno', rand(500, 8760));

   ice2lookup = createIce2Lookup();
end

function renamed = renameUsingIsmember(ice1, oldvars, newvars)
   % Rename using `ismember`.
   %
   % Historical note:
   %  - the function-based benchmark generally favored `ismember` slightly.
   %  - the old script-based benchmark once suggested `intersect` was nearly
   %    2x faster, but that result was treated cautiously because script scope
   %    and data reuse may have biased the measurement.
   renamed = renamevars(ice1, ...
      oldvars(ismember(oldvars, ice1.Properties.VariableNames)), ...
      newvars(ismember(oldvars, ice1.Properties.VariableNames)));
end

function renamed = renameUsingIntersect(ice1, oldvars, newvars)
   % Rename using `intersect`.
   %
   % Historical note:
   %  - this remained competitive and occasionally won on some runs, but the
   %    function-based benchmark usually made the gap look small enough that
   %    readability and consistency mattered more than the raw timing edge.
   [repvars, idx] = intersect(oldvars, ice1.Properties.VariableNames);
   renamed = renamevars(ice1, repvars, newvars(idx));
end

function ice2 = roundUsingSwitch(ice2)
   % Round using a `switch` over the available fields.
   %
   % Historical note:
   %  - this was one of the original generic round implementations.
   %  - this is also the pattern currently used in `icemodel.postprocess`
   %    inside `roundData(...)`.
   %  - the old real-data `time_roundData` comparison put the analogous
   %    switch-based postprocess rounder behind the direct-field variant but
   %    still close enough that the difference was mainly directional.
   fields = fieldnames(ice2);
   for mm = 1:numel(fields)
      thisfield = fields{mm};
      switch thisfield
         case {'f_ice', 'f_liq', 'k_vap', 'k_eff'}
            ice2.(thisfield) = round(ice2.(thisfield), 5);
         case {'Tice', 'h_melt', 'h_freeze'}
            ice2.(thisfield) = round(ice2.(thisfield), 3);
         case {'cp_sno', 'ro_sno'}
            ice2.(thisfield) = round(ice2.(thisfield), 1);
         case {'df_liq', 'df_lyr', 'Qsub', 'Sc', 'errT', 'errH'}
            ice2.(thisfield) = round(ice2.(thisfield), 8);
      end
   end
end

function ice2 = roundUsingIsmember(ice2, ice2lookup)
   % Round using `ismember` to select from the precision lookup.
   %
   % Historical note:
   %  - the old script-based benchmark often favored this rounding path.
   %  - the function-based benchmark usually kept it very close to the other
   %    generic lookup variants rather than showing a decisive win.
   lookup = ice2lookup(ismember(ice2lookup(:, 1), fieldnames(ice2)), :);
   for n = 1:size(lookup, 1)
      ice2.(lookup{n, 1}) = round(ice2.(lookup{n, 1}), lookup{n, 2});
   end
end

function ice2 = roundUsingIntersect(ice2, ice2lookup)
   % Round using `intersect` to select from the precision lookup.
   %
   % Historical note:
   %  - this often landed near or slightly ahead of the other generic lookup
   %    methods, but the spread was usually small enough to treat as a
   %    directional result rather than a decisive rule.
   [fields, idx] = intersect(ice2lookup(:, 1), fieldnames(ice2));
   for n = 1:numel(fields)
      ice2.(fields{n}) = round(ice2.(fields{n}), ice2lookup{idx(n), 2});
   end
end

function ice2 = roundUsingIsmemberPersistent(ice2)
   % Round using `ismember` with a persistent precision lookup.
   %
   % Historical note:
   %  - this variant came from the broader exploratory timing script rather
   %    than the original function-based benchmark.
   %  - it was preserved because the script-based comparison tested the
   %    persistent variants explicitly; in practice it stayed close to the
   %    non-persistent lookup methods.
   persistent ice2lookup
   if isempty(ice2lookup)
      ice2lookup = createIce2Lookup();
   end

   lookup = ice2lookup(ismember(ice2lookup(:, 1), fieldnames(ice2)), :);
   for n = 1:size(lookup, 1)
      ice2.(lookup{n, 1}) = round(ice2.(lookup{n, 1}), lookup{n, 2});
   end
end

function ice2 = roundUsingIntersectPersistent(ice2)
   % Round using `intersect` with a persistent precision lookup.
   %
   % Historical note:
   %  - this was also preserved from the exploratory timing script so the
   %    full set of compared generic rounders still exists in the formal
   %    benchmark.
   %  - like the other persistent variant, it usually stayed close to the
   %    non-persistent methods rather than clearly dominating them.
   persistent ice2lookup
   if isempty(ice2lookup)
      ice2lookup = createIce2Lookup();
   end

   [fields, idx] = intersect(ice2lookup(:, 1), fieldnames(ice2));
   for n = 1:numel(fields)
      ice2.(fields{n}) = round(ice2.(fields{n}), ice2lookup{idx(n), 2});
   end
end

function ice2lookup = createIce2Lookup()
%CREATEICE2LOOKUP Define the field-to-rounding-precision lookup table.
   ice2lookup = {
      'f_ice', 5; 'f_liq', 5; 'k_vap', 5; 'k_eff', 5; ...
      'Tice', 3; 'h_melt', 3; 'h_freeze', 3; ...
      'cp_sno', 1; 'ro_sno', 1; ...
      'df_liq', 8; 'df_lyr', 8; 'Qsub', 8; 'Sc', 8; 'errT', 8; ...
      'errH', 8};
end
