function [baseline_type, baseline_tag] = resolveBaselineSelector(baseline_tag)
   %RESOLVEBASELINESELECTOR Parse rolling vs release baseline selectors.
   %
   %  [baseline_type, baseline_tag] = ...
   %     icemodel.test.helpers.resolveBaselineSelector("rolling")
   arguments
      baseline_tag (1, :) string = "rolling"
   end

   % Empty strings and explicit "rolling" both select the mutable baseline.
   if isblanktext(baseline_tag) || strcmpi(baseline_tag, "rolling")
      baseline_type = "rolling";
      baseline_tag = "";
   else
      baseline_type = "release";
   end
end
