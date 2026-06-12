function [names, defaults] = metricRowSchema()
   %METRICROWSCHEMA Return canonical comparison-metric field names/defaults.
   %
   %  [names, defaults] = icemodel.verification.helpers.metricRowSchema()
   %
   % Outputs
   %  names      String column of metric table field names.
   %  defaults   Cell column of default values for one metric row.
   %
   % Role
   %  Operational helper that defines the comparecase metric table schema in
   %  one place so future metrics propagate consistently.

   % Keep the schema order aligned with comparecase output and summary CSVs.
   names = [ ...
      "experiment"
      "variable"
      "status"
      "n"
      "bias"
      "rmse"
      "correlation"
      "peak_target"
      "peak_candidate"
      "peak_error"
      "peak_time_error_hours"
      "snow_onset_time_error_hours"
      "melt_out_time_error_hours"];

   % Defaults represent a structurally valid row before any comparison occurs.
   defaults = { ...
      ""
      ""
      "ok"
      uint64(0)
      NaN
      NaN
      NaN
      NaN
      NaN
      NaN
      NaN
      NaN
      NaN};
end
