
% This is a custom script-based test, it does not conform to the performance
% testing framework

% Create a sample timetable
ice1 = timetable(datetime('now') + minutes(1:8760)', ...
   rand(8760, 1), rand(8760, 1), rand(8760, 1), rand(8760, 1), rand(8760, 1));
ice1.Properties.VariableNames = {'Qsi','Qsr','Qsn','Qli','Qle'};

% Define old and new variable names
oldvars = {'Qsi','Qsr','Qsn','Qli','Qle','Qln','Qh','Qe','Qc','Qn','Tsfc'};
newvars = {'swd','swu','swn','lwd','lwu','lwn','shf','lhf','chf','netr','tsfc'};

% Create a sample structure for rounding ice2 test
ice2 = struct('f_ice', rand(500, 8760), 'Tice', rand(500, 8760), ...
   'cp_sno', rand(500, 8760));

ice2lookup = {
   'f_ice', 5; 'f_liq', 5; 'k_vap', 5; 'k_eff', 5;
   'Tice', 3; 'h_melt', 3; 'h_freeze', 3;
   'cp_sno', 1; 'ro_sno', 1;
   'df_liq', 8; 'df_drn', 8; 'Qsub', 8; 'Sc', 8; 'errT', 8; 'errH', 8
   };

%% For initial testing, use tic toc

Rename Method 1: Using ismember
tmp = ice1;
tic
for m = 1:1000
tmp = renamevars(ice1, ...
   oldvars(ismember(oldvars, ice1.Properties.VariableNames)), ...
   newvars(ismember(oldvars, ice1.Properties.VariableNames)));
end
toc

% Rename Method 2: Using intersect
tic
for m = 1:1000
   [repvars, idx] = intersect(oldvars, ice1.Properties.VariableNames);
   tmp = renamevars(ice1, repvars, newvars(idx));
end
toc

%% Time rename

% Time both rename methods
f1 = @() rename1(oldvars, newvars, ice1);
f2 = @() rename2(oldvars, newvars, ice1);
time1 = timeit(f1);
time2 = timeit(f2);

fprintf('Time for rename method 1 (using ismember): %f seconds\n', time1);
fprintf('Time for rename method 2 (using intersect): %f seconds\n', time2);


% warmup / comparison using tic toc
tic
for m = 1:100
   lookup = ice2lookup(ismember(ice2lookup(:, 1), fieldnames(ice2)), :);
   for n = 1:size(lookup, 1)
      ice2.(lookup{n, 1}) = round(ice2.(lookup{n, 1}), lookup{n, 2});
   end
end
toc

tic
for m = 1:100
   [fields, idx] = intersect(ice2lookup(:, 1), fieldnames(ice2));
   for n = 1:numel(fields)
      ice2.(fields{n}) = round(ice2.(fields{n}), ice2lookup{idx(n), 2});
   end
end
toc


% Time rounding methods
f{1} = @() round1(ice2);
f{2} = @() round2(ice2, ice2lookup);
f{3} = @() round3(ice2, ice2lookup);
f{4} = @() round4(ice2);
f{5} = @() round5(ice2);

time = cell(numel(f), 1);
for n = 1:numel(f)
   time{n} = timeit(f{n});
end

methods = {'switch', 'ismember', 'intersect', 'ismember-persistent', 'intersect-persistent'};

for n = 1:numel(f)
   fprintf('Time for round method %d (using %s): %f seconds\n', n, methods{n}, time{n});
end

%% Local Functions

% Rename Method 1: Using ismember
function rename1(oldvars, newvars, ice1)
    ice1 = renamevars(ice1, ...
       oldvars(ismember(oldvars, ice1.Properties.VariableNames)), ...
       newvars(ismember(oldvars, ice1.Properties.VariableNames)));
end

% Rename Method 2: Using intersect
function rename2(oldvars, newvars, ice1)
    [repvars, idx] = intersect(oldvars, ice1.Properties.VariableNames);
    ice1 = renamevars(ice1, repvars, newvars(idx));
end

%%
% Round Method 1: switch case
function round1(ice2)
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
         case {'df_liq', 'df_drn', 'Qsub', 'Sc', 'errT', 'errH'}
            ice2.(thisfield) = round(ice2.(thisfield),8);
      end
   end
end

% Round Method 2: Using ismember
function round2(ice2, ice2lookup)
    lookup = ice2lookup(ismember(ice2lookup(:, 1), fieldnames(ice2)), :);
    for n = 1:size(lookup, 1)
        ice2.(lookup{n, 1}) = round(ice2.(lookup{n, 1}), lookup{n, 2});
    end
end

% Round Method 3: Using intersect
function round3(ice2, ice2lookup)
    [fields, idx] = intersect(ice2lookup(:, 1), fieldnames(ice2));
    for n = 1:numel(fields)
        ice2.(fields{n}) = round(ice2.(fields{n}), ice2lookup{idx(n), 2});
    end
end

% Round Method 4: Using ismember with persistent
function round4(ice2)
   persistent ice2lookup
   if isempty(ice2lookup)
      ice2lookup = {
         'f_ice', 5; 'f_liq', 5; 'k_vap', 5; 'k_eff', 5;
         'Tice', 3; 'h_melt', 3; 'h_freeze', 3;
         'cp_sno', 1; 'ro_sno', 1;
         'df_liq', 8; 'df_drn', 8; 'Qsub', 8; 'Sc', 8; 'errT', 8; 'errH', 8
         };
   end
   lookup = ice2lookup(ismember(ice2lookup(:, 1), fieldnames(ice2)), :);
   for n = 1:size(lookup, 1)
      ice2.(lookup{n, 1}) = round(ice2.(lookup{n, 1}), lookup{n, 2});
   end
end

% Round Method 5: Using intersect with persistent
function round5(ice2)
   persistent ice2lookup
   if isempty(ice2lookup)
      ice2lookup = {
         'f_ice', 5; 'f_liq', 5; 'k_vap', 5; 'k_eff', 5;
         'Tice', 3; 'h_melt', 3; 'h_freeze', 3;
         'cp_sno', 1; 'ro_sno', 1;
         'df_liq', 8; 'df_drn', 8; 'Qsub', 8; 'Sc', 8; 'errT', 8; 'errH', 8
         };
   end
    [fields, idx] = intersect(ice2lookup(:, 1), fieldnames(ice2));
    for n = 1:numel(fields)
        ice2.(fields{n}) = round(ice2.(fields{n}), ice2lookup{idx(n), 2});
    end
end

