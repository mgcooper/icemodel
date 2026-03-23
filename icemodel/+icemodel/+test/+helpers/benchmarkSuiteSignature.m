function [signature, files] = benchmarkSuiteSignature()
   %BENCHMARKSUITESIGNATURE Hash the managed core benchmark suite.
   %
   %  signature = icemodel.test.helpers.benchmarkSuiteSignature()
   %  [signature, files] = icemodel.test.helpers.benchmarkSuiteSignature()
   %
   % This signature is used to decide whether a benchmark baseline saved
   % inside a perf baseline file is still comparable to the currently
   % checked-out benchmark suite.

   % Collect the public runner plus the default top-level benchmark files.
   testdir = icemodel.getpath('test');
   benchdir = fullfile(testdir, 'benchmarks');
   listing = dir(fullfile(benchdir, '*.m'));
   files = string(sort({listing.name}.'));
   files = fullfile(benchdir, files);
   files = vertcat(fullfile(testdir, 'run_benchmark_suite.m'), files);

   % Hash the file names and contents so changed benchmark bodies or runner
   % selection logic invalidate stale embedded benchmark baselines.
   md = java.security.MessageDigest.getInstance('MD5');
   for i = 1:numel(files)
      pathname = char(files(i));
      md.update(uint8(pathname));
      md.update(uint8(10));
      md.update(uint8(fileread(pathname)));
      md.update(uint8(10));
   end

   % Format the digest as a stable lower-case hex string.
   digest = typecast(md.digest(), 'uint8');
   signature = string(lower(reshape(dec2hex(digest, 2).', 1, [])));
end
