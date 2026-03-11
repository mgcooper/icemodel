function S = runPerfCase(experiment, suite, c)
   %RUNPERFCASE Run one formal performance case and normalize the result.
   %
   %  S = test.helpers.runPerfCase(experiment, suite, c)
   %
   % Output fields:
   %  samples, activity, sample_times, activity_times, valid, n_warmups

   configurePerfCaseEnv(c);
   perf_result = run(experiment, suite);
   perf_result = perf_result(1);

   samples = resultTable(perf_result, 'Samples');
   activity = resultTable(perf_result, 'TestActivity');
   tcol = measuredTimeColumn(samples);
   acol = measuredTimeColumn(activity);
   if isempty(tcol)
      error('unable to find measured time column in perf result Samples table')
   end

   S = struct();
   S.samples = samples;
   S.activity = activity;
   S.sample_times = toSeconds(samples.(tcol));
   S.valid = logical(perf_result.Valid);
   S.n_warmups = countWarmups(activity);

   if isempty(acol)
      S.activity_times = nan(height(activity), 1);
   else
      S.activity_times = toSeconds(activity.(acol));
   end
end

function configurePerfCaseEnv(c)
   setenv('ICEMODEL_TEST_SMBMODEL', char(c.smbmodel));
   setenv('ICEMODEL_TEST_SITENAME', char(c.sitename));
   setenv('ICEMODEL_TEST_FORCINGS', char(c.forcings));
   setenv('ICEMODEL_TEST_USERDATA', char(c.userdata));
   setenv('ICEMODEL_TEST_USERVARS', char(c.uservars));
   setenv('ICEMODEL_TEST_SIMYEAR', int2str(c.simyear));
   setenv('ICEMODEL_TEST_SOLVER', int2str(c.solver));
end

function T = resultTable(result, propname)
   T = result.(propname);
   if isstruct(T)
      T = struct2table(T);
   end
end

function name = measuredTimeColumn(T)
   name = '';
   if isempty(T)
      return
   end
   names = string(T.Properties.VariableNames);
   cand = ["MeasuredTime", "MeasuredValue", "WallTime"];
   idx = find(ismember(cand, names), 1, 'first');
   if ~isempty(idx)
      name = char(cand(idx));
   end
end

function x = toSeconds(x)
   if isduration(x)
      x = seconds(x);
   end
end

function n = countWarmups(activity)
   n = 0;
   if isempty(activity) || ...
         ~ismember('Objective', activity.Properties.VariableNames)
      return
   end
   n = sum(strcmpi(string(activity.Objective), 'Warmup'));
end
