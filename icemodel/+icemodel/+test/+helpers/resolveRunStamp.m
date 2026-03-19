function [run_date, run_id, run_name] = resolveRunStamp(run_name)
   %RESOLVERUNSTAMP Resolve shared batch run identifiers for test artifacts.
   %
   %  [run_date, run_id, run_name] = ...
   %     icemodel.test.helpers.resolveRunStamp(run_name)
   arguments
      run_name (1, :) string = ""
   end

   % Generate a local timestamp when the caller did not request one.
   if isblanktext(run_name)
      t = datetime('now', 'TimeZone', 'local');
      run_date = string(datetime(t, 'Format', 'yyyyMMdd'));
      run_id = string(datetime(t, 'Format', 'HHmmss'));
      run_name = run_date + "-" + run_id;
      return
   end

   % Accept only the canonical yyyymmdd-HHMMSS run naming format.
   parts = split(run_name, "-");
   if numel(parts) ~= 2
      error('run_name must have format yyyymmdd-HHMMSS')
   end
   run_date = parts(1);
   run_id = parts(2);
end
