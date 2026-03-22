clear functions

% Toggle the formal perf-suite call at the bottom of the script.
run_perf_suite_flag = false;

% Build one resolved smoke-style options struct for the direct timing sweep.
opts0 = icemodel.setopts();
opts0 = icemodel.resetopts(opts0, ...
   'smbmodel', 'icemodel', ...
   'sitename', 'kanm', ...
   'solver', 2, ...
   'simyears', 2016, ...
   'n_spinup_years', 0);
opts0 = icemodel.configureRun(opts0);

% Time each spectral variant directly around the main model call.
variants = ["inlined"; "functions"; "lookup"];
times = nan(size(variants));

for i = 1:numel(variants)
   opts = opts0;
   opts.test_spectral_variant = variants(i);
   opts.testname = "manual_" + variants(i);

   icemodel(opts); % warmup

   t0 = tic;
   icemodel(opts);
   times(i) = toc(t0);

   fprintf('%s: %.3f s\n', variants(i), times(i));
end

table(variants, times)

%%

% Optionally compare the same accepted lookup path through the formal suite.
if run_perf_suite_flag

   clear functions
   r = run_perf_suite( ...
      tier="smoke", ...
      smbmodel="icemodel", ...
      solver=2, ...
      baseline="rolling", ...
      n_runs=1, ...
      include_benchmarks=false, ...
      spectral_variant="lookup");
end
