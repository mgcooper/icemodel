function [run_date, run_id, run_name] = resolveRunStamp(run_name)
%RESOLVERUNSTAMP Resolve shared batch run identifiers for test artifacts.
   arguments
      run_name = string.empty()
   end

   run_name = string(run_name);

   if isempty(run_name) || all(strlength(run_name) == 0)
      t = datetime('now', 'TimeZone', 'local');
      run_date = string(datestr(t, 'yyyymmdd'));
      run_id = string(datestr(t, 'HHMMSS'));
      run_name = run_date + "-" + run_id;
      return
   end

   parts = split(run_name, "-");
   if numel(parts) ~= 2
      error('run_name must have format yyyymmdd-HHMMSS')
   end
   run_date = parts(1);
   run_id = parts(2);
end
