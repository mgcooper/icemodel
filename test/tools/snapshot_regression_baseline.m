function RegressionBaseline = snapshot_regression_baseline(kwargs)
   %SNAPSHOT_REGRESSION_BASELINE Save a release snapshot from the rolling regression baseline.
   %
   %  RegressionBaseline = snapshot_regression_baseline(baseline_tag="v1.1")
   %
   % Use this when the current rolling regression baseline should be frozen
   % as a named release baseline. This does not rerun the model; it copies
   % the current rolling baseline into a versioned release file.

   arguments (Input)
      kwargs.baseline_tag (1, :) string
      kwargs.smbmodel (1, :) string {mustBeMember(kwargs.smbmodel, ...
         ["all", "icemodel", "skinmodel"])} = "all"
      kwargs.overwrite (1, 1) logical = false
      kwargs.output_file string = string.empty()
   end
   [baseline_tag, smbmodel, overwrite, output_file] = deal( ...
      kwargs.baseline_tag, kwargs.smbmodel, kwargs.overwrite, ...
      kwargs.output_file);

   % Copy the current rolling regression baseline into a versioned release file.
   RegressionBaseline = test.helpers.snapshotBaseline( ...
      "regression", baseline_tag, smbmodel, overwrite, output_file);
end
