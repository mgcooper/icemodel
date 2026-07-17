function name = findEnclosingWindowFile(directory, prefix, suffix, qstart, qend, kwargs)
   %FINDENCLOSINGWINDOWFILE Name of a staged window file bracketing a query span.
   %
   %  name = icemodel.forcing.helpers.findEnclosingWindowFile(DIRECTORY, ...
   %     PREFIX, SUFFIX, QSTART, QEND)
   %
   % Returns the file name PREFIX_<YYYYMMDD>_<YYYYMMDD>SUFFIX in DIRECTORY whose
   % encoded period [file-start 00:00, file-end 23:59:59] CONTAINS the query span
   % [QSTART, QEND], or "" when none matches. When several files enclose the
   % query, selection is deterministic: widest window first, then latest end,
   % then lexical file name. This is the single source of the
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

   arguments
      directory (1, 1) string
      prefix (1, 1) string
      suffix (1, 1) string
      qstart (1, 1) datetime
      qend (1, 1) datetime
      kwargs.accept_candidate = []
   end

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
   candidates = strings(numel(d), 1);
   starts = NaT(numel(d), 1, 'TimeZone', tz);
   ends = NaT(numel(d), 1, 'TimeZone', tz);
   n_candidates = 0;
   for n = 1:numel(d)
      tok = regexp(d(n).name, pat, 'tokens', 'once');
      if isempty(tok)
         continue
      end
      file_start = datetime(tok{1}, 'InputFormat', 'yyyyMMdd', 'TimeZone', tz);
      file_end = datetime(tok{2}, 'InputFormat', 'yyyyMMdd', 'TimeZone', tz) ...
         + hours(23) + minutes(59) + seconds(59);
      if file_start <= qstart && file_end >= qend
         candidate_path = fullfile(directory, string(d(n).name));
         if ~isempty(kwargs.accept_candidate) ...
               && ~kwargs.accept_candidate(candidate_path)
            continue
         end
         n_candidates = n_candidates + 1;
         candidates(n_candidates) = string(d(n).name);
         starts(n_candidates) = file_start;
         ends(n_candidates) = file_end;
      end
   end
   if n_candidates == 0
      return
   end

   % Select by explicit sortable keys instead of filesystem directory order.
   candidates = candidates(1:n_candidates);
   starts = starts(1:n_candidates);
   ends = ends(1:n_candidates);
   durations = seconds(ends - starts);
   end_key = year(ends) .* 10000 + month(ends) .* 100 + day(ends);
   rank = table(-durations, -end_key, candidates, ...
      'VariableNames', {'neg_duration', 'neg_end', 'name'});
   [~, order] = sortrows(rank, {'neg_duration', 'neg_end', 'name'});
   name = candidates(order(1));
end
