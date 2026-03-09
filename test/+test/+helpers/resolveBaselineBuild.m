function [baseline_type, baseline_tag, output_file] = resolveBaselineBuild( ...
      kind, baseline, baseline_tag, smbmodel, output_file, simyear)
%RESOLVEBASELINEBUILD Resolve baseline type/tag and default output file.
%
%  [baseline_type, baseline_tag, output_file] = ...
%     test.helpers.resolveBaselineBuild("perf", baseline, baseline_tag, smbmodel, output_file, simyear)

   arguments
      kind (1, :) string {mustBeMember(kind, ["perf", "regression"])}
      baseline string = string.empty()
      baseline_tag string = string.empty()
      smbmodel string = "all"
      output_file string = string.empty()
      simyear double = NaN
   end

   if isblanktext(baseline)
      baseline = baseline_tag;
   end
   if isblanktext(baseline)
      error('baseline or baseline_tag is required')
   end

   if baseline == "rolling"
      baseline_type = "rolling";
      baseline_tag = "";
   else
      baseline_type = "release";
      baseline_tag = baseline;
   end

   if isblanktext(output_file)
      output_file = test.helpers.defaultBaselinePath( ...
         kind, baseline_type, baseline_tag, smbmodel, simyear);
   end
end
