function [t1, t2, t3, t4] = time_roundData(ice1, ice2)
% based on these tests they're all pretty similar, round2 has a slight edge,
% then round3 (persistent), then round4. Clearly eliminate round1 in favor of
% round2 

% For these tests, need to run icemodel with the roudn step turned off in
% POSTPROC, which is annoying but its the actual data, see test_renameRound for
% synthetic testing.

ice1test = ice1;
ice2test = ice2;

% warm up on the og data
[ice1, ice2] = round1(ice1, ice2);
[ice1, ice2] = round2(ice1, ice2, opts);
[ice1, ice2] = round3(ice1, ice2);
[ice1, ice2] = round4(ice1, ice2);

f1 = @() round1(ice1test, ice2test);
f2 = @() round2(ice1test, ice2test, opts);
f3 = @() round3(ice1test, ice2test);
f4 = @() round4(ice1test, ice2test);

t1 = timeit(f1, 2)
t2 = timeit(f2, 2)
t3 = timeit(f3, 2)
t4 = timeit(f4, 2)


%% round1
function [ice1, ice2] = round1(ice1, ice2)

% Round ice1 data to 5 digits (all fields, so there are no if-else checks).
fields = ice1.Properties.VariableNames;
for mm = 1:numel(fields)
   thisfield = fields{mm};
   if ismember(thisfield,fields)
      ice1.(thisfield) = round(ice1.(thisfield),5);
   end
end

% Round the ice2 data
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

%% round2
function [ice1, ice2] = round2(ice1, ice2, opts)

% round. this should be faster than the two options below
ice1.Tsfc = round(ice1.Tsfc,5);
ice1.melt = round(ice1.melt,5);
ice1.runoff = round(ice1.runoff,5);
ice1.freeze = round(ice1.freeze,5);
try
   ice1.drain = round(ice1.drain,5);
catch
end

ice2.Tice = round(ice2.Tice,3);
if strcmp('icemodel', opts.simmodel)
   ice2.f_ice = round(ice2.f_ice,5);
   ice2.f_liq = round(ice2.f_liq,5);
   ice2.df_liq = round(ice2.df_liq,8);
end
try
   ice2.df_drn = round(ice2.df_drn,8);
catch
end

%% round3
function [ice1, ice2] = round3(ice1, ice2)

% Round the ice1 data to five digits
ice1{:, ice1.Properties.VariableNames} = ...
   round(ice1{:, ice1.Properties.VariableNames}, 5);

% Round the ice2 data to variable precisioin
persistent ice2lookup
if isempty(ice2lookup)
   ice2lookup = {
       'f_ice', 5; 'f_liq', 5; 'k_vap', 5; 'k_eff', 5;
       'Tice', 3; 'h_melt', 3; 'h_freeze', 3;
       'cp_sno', 1; 'ro_sno', 1;
       'df_liq', 8; 'df_drn', 8; 'Qsub', 8; 'Sc', 8; 'errT', 8; 'errH', 8
   };
end

% Apply the rounding to ice2 according to the defined precision
lookup = ice2lookup(ismember(ice2lookup(:, 1), fieldnames(ice2)), :);
for n = 1:size(lookup, 1)
    ice2.(lookup{n, 1}) = round(ice2.(lookup{n, 1}), lookup{n, 2});
end

%% round4
function [ice1, ice2] = round4(ice1, ice2)

% Round the ice1 data to five digits
ice1{:, ice1.Properties.VariableNames} = ...
   round(ice1{:, ice1.Properties.VariableNames}, 5);

% Round the ice2 data
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


