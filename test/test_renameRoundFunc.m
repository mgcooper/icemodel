function tests = test_renameRoundFunc
   tests = functiontests(localfunctions);

   %#ok<*INUSD>
   %#ok<*NASGU>
end

% function setup(testCase)
%
%    % Save the test data
%    [ice1, oldvars, newvars] = createIce1();
%
%    testCase.TestData = generateTestData(funcname);
%
%    % example test data
%    testCase.TestData.S = struct('field1', 'value1', 'field2', 'value2');
%
%    % see various built in "Fixture" functions to perform common setup and
%    % teardown tasks. You can also write your own, for example to formalize
%    % generateTestData
% end

% test rename using ismember
function testRenameIsmember(testCase)
   [ice1, oldvars, newvars] = createIce1();
   ice1 = renamevars(ice1, ...
      oldvars(ismember(oldvars, ice1.Properties.VariableNames)), ...
      newvars(ismember(oldvars, ice1.Properties.VariableNames)));
end

% test rename using intersect
function testRenameIntersect(testCase)

   [ice1, oldvars, newvars] = createIce1();

   [repvars, idx] = intersect(oldvars, ice1.Properties.VariableNames);
   ice1 = renamevars(ice1, repvars, newvars(idx));
end

% test round using switch case
function testRoundSwitch(testCase)

   [ice2, ~] = createIce2();

   fields = fieldnames(ice2);
   for mm = 1:numel(fields)
      thisfield = fields{mm};
      switch thisfield
         case {'f_ice', 'f_liq', 'k_vap', 'k_eff'}
            ice2.(thisfield) = round(ice2.(thisfield),5);
         case {'Tice', 'h_melt', 'h_freeze'}
            ice2.(thisfield) = round(ice2.(thisfield),3);
         case {'cp_sno', 'ro_sno'}
            ice2.(thisfield) = round(ice2.(thisfield),1);
         case {'df_liq', 'df_lyr', 'Qsub', 'Sc', 'errT', 'errH'}
            ice2.(thisfield) = round(ice2.(thisfield),8);
      end
   end
end

% test round using ismember
function testRoundIsmember(testCase)

   [ice2, ice2lookup] = createIce2();

   lookup = ice2lookup(ismember(ice2lookup(:, 1), fieldnames(ice2)), :);
   for n = 1:size(lookup, 1)
      ice2.(lookup{n, 1}) = round(ice2.(lookup{n, 1}), lookup{n, 2});
   end
end

% test round using intersect
function testRoundIntersect(testCase)

   [ice2, ice2lookup] = createIce2();

   [fields, idx] = intersect(ice2lookup(:, 1), fieldnames(ice2));
   for n = 1:numel(fields)
      ice2.(fields{n}) = round(ice2.(fields{n}), ice2lookup{idx(n), 2});
   end
end

% create ice1 test data
function [ice1, oldvars, newvars] = createIce1()
   % Create a sample timetable
   ice1 = timetable(datetime('now') + minutes(1:8760)', ...
      rand(8760, 1), rand(8760, 1), rand(8760, 1), rand(8760, 1), rand(8760, 1));
   ice1.Properties.VariableNames = {'Qsi','Qsr','Qsn','Qli','Qle'};

   % Define old and new variable names
   oldvars = {'Qsi','Qsr','Qsn','Qli','Qle','Qln','Qh','Qe','Qc','Qn','Tsfc'};
   newvars = {'swd','swu','swn','lwd','lwu','lwn','shf','lhf','chf','netr','tsfc'};
end

% create ice2 test data
function [ice2, ice2lookup] = createIce2()
   % Create a sample structure for rounding ice2 test
   ice2 = struct('f_ice', rand(500, 8760), 'Tice', rand(500, 8760), ...
      'cp_sno', rand(500, 8760));

   ice2lookup = {
      'f_ice', 5; 'f_liq', 5; 'k_vap', 5; 'k_eff', 5;
      'Tice', 3; 'h_melt', 3; 'h_freeze', 3;
      'cp_sno', 1; 'ro_sno', 1;
      'df_liq', 8; 'df_lyr', 8; 'Qsub', 8; 'Sc', 8; 'errT', 8; 'errH', 8
      };
end
