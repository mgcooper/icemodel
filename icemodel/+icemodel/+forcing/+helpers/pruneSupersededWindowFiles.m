function removed = pruneSupersededWindowFiles(new_file, prefix, suffix, kwargs)
   %PRUNESUPERSEDEDWINDOWFILES Remove shorter windows contained by a new file.
   %
   %  removed = icemodel.forcing.helpers.pruneSupersededWindowFiles( ...
   %     new_file, prefix, suffix)
   %
   % Call only after NEW_FILE was written successfully. Files in the same
   % directory and naming/cadence class are removed when their encoded date
   % window is strictly contained by the new window. Overlapping, enclosing,
   % unrelated, and equal-window files remain untouched. A warning lists every
   % removal before deletion so wider refresh side effects are explicit.

   arguments
      new_file (1, 1) string
      prefix (1, 1) string
      suffix (1, 1) string
      kwargs.accept_candidate = []
   end

   removed = strings(0, 1);
   [directory, name, extension] = fileparts(new_file);
   new_name = string(name) + string(extension);
   [ok, new_start, new_end] = parseWindowName(new_name, prefix, suffix);
   if ~ok || ~isfolder(directory)
      return
   end

   % Collect strictly contained candidates in lexical order for deterministic
   % warning and deletion behavior.
   listing = dir(fullfile(directory, char(prefix + "_*_*" + suffix)));
   names = sort(string({listing.name}));
   candidates = strings(numel(names), 1);
   n_candidates = 0;
   for k = 1:numel(names)
      if names(k) == new_name
         continue
      end
      [valid, candidate_start, candidate_end] = ...
         parseWindowName(names(k), prefix, suffix);
      contained = valid && candidate_start >= new_start ...
         && candidate_end <= new_end ...
         && (candidate_start > new_start || candidate_end < new_end);
      candidate_path = fullfile(directory, names(k));
      if contained && ~isempty(kwargs.accept_candidate)
         contained = kwargs.accept_candidate(candidate_path);
      end
      if contained
         n_candidates = n_candidates + 1;
         candidates(n_candidates) = candidate_path;
      end
   end
   removed = candidates(1:n_candidates);
   if isempty(removed)
      return
   end

   % Warn before the irreversible cleanup, then delete exactly the listed files.
   warning('icemodel:forcing:pruneSupersededWindowFiles:removed', ...
      ['Wider artifact %s supersedes contained shorter artifact(s): %s. ' ...
      'Removing the shorter file(s).'], new_file, strjoin(removed, ', '));
   for filename = reshape(removed, 1, [])
      delete(filename);
   end
end

function [ok, window_start, window_end] = ...
      parseWindowName(name, prefix, suffix)
   %PARSEWINDOWNAME Decode one PREFIX_YYYYMMDD_YYYYMMDD_SUFFIX file name.
   pattern = ['^' regexptranslate('escape', char(prefix)) ...
      '_(\d{8})_(\d{8})' regexptranslate('escape', char(suffix)) '$'];
   token = regexp(char(name), pattern, 'tokens', 'once');
   ok = ~isempty(token);
   window_start = NaT('TimeZone', 'UTC');
   window_end = NaT('TimeZone', 'UTC');
   if ~ok
      return
   end
   window_start = datetime(token{1}, 'InputFormat', 'yyyyMMdd', ...
      'TimeZone', 'UTC');
   window_end = datetime(token{2}, 'InputFormat', 'yyyyMMdd', ...
      'TimeZone', 'UTC') + days(1) - seconds(1);
end
