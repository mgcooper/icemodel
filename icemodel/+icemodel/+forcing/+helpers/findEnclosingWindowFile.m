function name = findEnclosingWindowFile(directory, prefix, suffix, qstart, qend)
   %FINDENCLOSINGWINDOWFILE Name of a staged window file bracketing a query span.
   %
   %  name = icemodel.forcing.helpers.findEnclosingWindowFile(DIRECTORY, ...
   %     PREFIX, SUFFIX, QSTART, QEND)
   %
   % Returns the file name PREFIX_<YYYYMMDD>_<YYYYMMDD>SUFFIX in DIRECTORY whose
   % encoded period [file-start 00:00, file-end 23:59:59] CONTAINS the query span
   % [QSTART, QEND], or "" when none matches. This is the single source of the
   % "does a staged full-period file already cover this window" logic shared by
   % met-file resolution (icemodel.createMetFileNames) and met-swap userdata
   % resolution (icemodel.loadmet).
   %
   % Inputs
   %  directory - folder to search
   %  prefix    - file name prefix up to (not including) the date pair, e.g.
   %              "met_kanm_promice" or "kanm_racmo"
   %  suffix    - file name suffix after the date pair, e.g. "_1hr.mat" or ".mat"
   %  qstart    - datetime, start of the requested span (its TimeZone is applied
   %              to the parsed file dates so the comparison is valid)
   %  qend      - datetime, end of the requested span
   %
   % Output
   %  name      - the matching file name (string), or "" if none brackets

   name = "";
   if ~isfolder(directory)
      return
   end
   prefix = char(prefix);
   suffix = char(suffix);
   d = dir(fullfile(char(directory), [prefix '_*_*' suffix]));
   pat = ['^' regexptranslate('escape', prefix) ...
      '_(\d{8})_(\d{8})' regexptranslate('escape', suffix) '$'];
   tz = qstart.TimeZone;
   for n = 1:numel(d)
      tok = regexp(d(n).name, pat, 'tokens', 'once');
      if isempty(tok)
         continue
      end
      file_start = datetime(tok{1}, 'InputFormat', 'yyyyMMdd', 'TimeZone', tz);
      file_end = datetime(tok{2}, 'InputFormat', 'yyyyMMdd', 'TimeZone', tz) ...
         + hours(23) + minutes(59) + seconds(59);
      if file_start <= qstart && file_end >= qend
         name = string(d(n).name);
         return
      end
   end
end
