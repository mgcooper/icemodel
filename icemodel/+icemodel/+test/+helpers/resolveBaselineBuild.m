function [baseline_type, baseline_tag, output_file] = resolveBaselineBuild( ...
      kind, baseline, baseline_tag, smbmodel, output_file, simyear)
   %RESOLVEBASELINEBUILD Resolve baseline type/tag and default output file.
   %
   %  [baseline_type, baseline_tag, output_file] = ...
   %     icemodel.test.helpers.resolveBaselineBuild("perf", baseline, baseline_tag, smbmodel, output_file, simyear)

   arguments
      kind (1, :) string {mustBeMember(kind, ["perf", "regression"])}
      baseline string = string.empty()
      baseline_tag string = string.empty()
      smbmodel string = "all"
      output_file string = string.empty()
      simyear double = NaN
   end

   % Default the public build API to the managed rolling baseline.
   if isblanktext(baseline) && isblanktext(baseline_tag)
      baseline = "rolling";
   elseif isblanktext(baseline)
      baseline = baseline_tag;
   elseif baseline == "rolling" && ~isblanktext(baseline_tag)
      % Let BASELINE_TAG upgrade the default rolling selector into a
      % concrete release target without requiring the caller to blank
      % BASELINE first.
      baseline = baseline_tag;
   end

   % Resolve the selector into the canonical baseline type/tag pair.
   if baseline == "rolling"
      baseline_type = "rolling";
      baseline_tag = "";
   else
      baseline_type = "release";
      baseline_tag = baseline;
   end

   % Fill the default managed output path when the caller did not request
   % a custom file.
   if isblanktext(output_file)
      output_file = icemodel.test.helpers.defaultBaselinePath( ...
         kind, baseline_type, baseline_tag, smbmodel, simyear);
   end
end
