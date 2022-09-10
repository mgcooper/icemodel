
% TLDR: v7.3 is the winner. next need to check time to process and convert
% to timetable versus 

% TLDR: v6 saves about 10x faster but is also 10x larger. v4 is even faster
% than v6 but same size. Therefore if we really want speed and only save
% arrays, we can use v4, but we need tons of space. This would make sense
% for running lots of simulations where the possibility of re-running is
% high, and then doing compression after we're finished. 

% file type          MB          tt (s)   struct (s)  vec (s)
% v4:                141         nan      nan         0.39
% v6:                141         0.45     0.52        0.49
% v7:                13.5        2.26     2.21        2.18
% v7.3:              11          1.07     1.04        1.06
% v7, no comp:       141         1.01     1.01        1.05
% v7.3, no comp:     141         0.26     0.33        0.38

% in terms of speed across the different POSTPROCS:

% POSTPROC2          2.66
% POSTPROC_struct    1.22
% POSTPROC_vec       1.15

% but this isn't entirely fair bcuase POSTPROC does the unneccessary
% runoff/h_melt/h_freeze stuff because I cannot edit it until the runs are
% finished, so i need to repeat this test



% below here is basically summarized by the table above

% Across formats, file sizes are identical i.e. v6 struct is same as v6 tt
% and v6 array (as long as the hourly arrays are saved, not the raw 15 min)

% Across versions, file sizes are about 10x larger for v4/v6 (141 vs 13 MB)
% so we can rule them out based on size alone, since there also is not a
% huge increase in speed (see next)

% Across formats, v6 is about 4x faster than v7 (1/2 sec vs 2.5 sec) and
% about 2x faster than v7.3 (1/2 sec vs 1.1 sec)

% Across versions, v6 is about 4x faster than v7 (1/2 sec vs 2.5 sec) and
% about 2x faster than v7.3 (1/2 sec vs 1.1 sec)

% in terms of speed, there is only a slight advantage to using v6 relative
% to v7.3, but both v6 and v7.3 are about 2x as fast as v7. Specifically
% v7.3 is 1.45 sec/year/point whereas v7 is 2.85 sec/year/point. For 10
% years and 1487 points, that could save 6 hours. 

% THEREFORE, it seems clear that the best approach is v7.3, with whatever
% format is most convenient for future steps i.e. timetable, struct, or
% arrays


% saving the raw output is not a good option b/c it is 15 min and we need
% the high precision to compute derived vars, so we want to do the post
% processing (as limited as possible) before saving

% saving as arrays rather than timetables might be a good idea, esp for the
% gridded sector runs, since we need to compute basin-scale values later

% 
clean
load('test/tt/times'); times_tt = times; clear times;
load('test/struct/times'); times_struct = times; clear times;
load('test/vec/times'); times_vec = times; clear times;

% % v4 tt and struct nothing was saved - structs and tt's are not supported
% % for tt and struct, the 'out' row can be discarded b/c we know for sure we
% % don't want to save the raw 15 min, too much data. Finally, to compare the
% % tt and struct to the arrays, we should add ice1+ice2 time. therefore we
% % make the following changes:
% times_tt.v4(1) = nan;
% times_tt = times_tt(1:2,:);
% times_struct.v4(1) = nan;
% times_struct = times_struct(1:2,:);
% 
% % now we can compare them, but we should add
% test_tt = sum(table2array(times_tt));
% test_struct = sum(table2array(times_struct));
% test_vec = table2array(times_vec);
% 
% times = [test_tt;test_struct;test_vec];
% times = array2table(times,'VariableNames',{'v4','v6','v7','v73'}, ...
%    'RowNames',{'timetable','struct','arrays'})

