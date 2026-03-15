function [baseline_type, baseline_tag] = resolveBaselineSelector(baseline_tag)
%RESOLVEBASELINESELECTOR Parse rolling vs release perf baseline selector.
   arguments
      baseline_tag = "rolling"
   end

   baseline_tag = string(baseline_tag);
   if isempty(baseline_tag) || ...
         (isstring(baseline_tag) && all(strlength(baseline_tag) == 0)) || ...
         lower(string(baseline_tag)) == "rolling"
      baseline_type = "rolling";
      baseline_tag = "";
   else
      baseline_type = "release";
   end
end
