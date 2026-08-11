function [series, filename] = loadWidestTimetable(hits)
   %LOADWIDESTTIMETABLE Load the staged timetable with the widest time axis.
   %
   %  [series, filename] = ...
   %     icemodel.forcing.reconstruct.loadWidestTimetable(dir(...))
   %
   % MAT-file size is not a measure of coverage: metadata and unrelated
   % payloads can make a narrow window larger. The function selects the
   % greatest saved time span and breaks ties on row count. The function skips
   % files without a nonempty timetable.
   % filename is the selected artifact path, or "" when no file is usable.

   arguments
      hits (1, :) struct
   end

   series = [];
   filename = "";
   best_span = -Inf;
   best_height = -Inf;
   for k = 1:numel(hits)
      S = load(fullfile(hits(k).folder, hits(k).name));
      vars = fieldnames(S);
      for v = 1:numel(vars)
         candidate = S.(vars{v});
         if ~istimetable(candidate) || isempty(candidate)
            continue
         end

         % A single-row timetable has zero span and remains a valid fallback.
         span = 0;
         if height(candidate) > 1
            span = seconds(candidate.Properties.RowTimes(end) ...
               - candidate.Properties.RowTimes(1));
         end
         if span > best_span || ...
               (span == best_span && height(candidate) > best_height)
            series = candidate;
            filename = string(fullfile(hits(k).folder, hits(k).name));
            best_span = span;
            best_height = height(candidate);
         end
      end
   end
end
