function PerfBaseline = snapshot_perf_baseline(kwargs)
   %SNAPSHOT_PERF_BASELINE Save release snapshot from the rolling perf baseline.
   %
   %  PerfBaseline = snapshot_perf_baseline(baseline_tag="v1.1")
   %  PerfBaseline = snapshot_perf_baseline(baseline_tag="v1.1", simyear=2016)
   %
   % Use this when the current rolling perf baseline should be frozen as a
   % named release baseline. This does not rerun the model; it copies the
   % current rolling baseline into a versioned release file.

   arguments (Input)
      kwargs.baseline_tag (1, :) string
      kwargs.smbmodel (1, :) string {mustBeMember(kwargs.smbmodel, ...
         ["all", "icemodel", "skinmodel"])} = "all"
      kwargs.simyear (1, 1) double {mustBeInteger, mustBePositive} = 2016
      kwargs.overwrite (1, 1) logical = false
      kwargs.output_file string = string.empty()
   end
   [baseline_tag, smbmodel, simyear, overwrite, output_file] = deal( ...
      kwargs.baseline_tag, kwargs.smbmodel, kwargs.simyear, ...
      kwargs.overwrite, kwargs.output_file);

   % Copy the current rolling perf baseline into a versioned release file.
   PerfBaseline = test.helpers.snapshotBaseline( ...
      "perf", baseline_tag, smbmodel, overwrite, output_file, simyear);
end
