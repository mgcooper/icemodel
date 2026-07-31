function version = readCffVersion(cffpath)
   %READCFFVERSION Read the top-level software version from a CFF file.
   %
   %  VERSION = ICEMODEL.INTERNAL.READCFFVERSION(CFFPATH) returns the
   %  top-level version scalar as a character row vector. The scalar may be
   %  unquoted, single quoted, or double quoted.
   %
   % See also: ICEMODEL.INTERNAL.VERSION

   % Accept the character paths used by the legacy runtime as well as scalar
   % strings without requiring syntax newer than the supported MATLAB floor.
   if isstring(cffpath)
      if ~isscalar(cffpath)
         error('icemodel:internal:readCffVersion:InvalidPath', ...
            'CFF path must be a text scalar.');
      end
      cffpath = char(cffpath);
   elseif ~(ischar(cffpath) && isrow(cffpath))
      error('icemodel:internal:readCffVersion:InvalidPath', ...
         'CFF path must be a text scalar.');
   end

   % Fail before reading so a missing release-archive citation file has a
   % stable, actionable error rather than a generic file I/O failure.
   if exist(cffpath, 'file') ~= 2
      error('icemodel:internal:readCffVersion:FileNotFound', ...
         'CFF file not found: %s', cffpath);
   end

   % Match only a top-level key because nested version fields must not become
   % the persisted IceModel version by accident. A YAML mapping separator
   % requires whitespace after the colon unless the value is empty.
   contents = fileread(cffpath);
   lines = regexp(contents, '\r\n|\n|\r', 'split');
   matches = lines(~cellfun('isempty', regexp(lines, ...
      '^version[ \t]*:(?:[ \t]+.*|[ \t]*)$', 'once')));
   if isempty(matches)
      error('icemodel:internal:readCffVersion:MissingVersion', ...
         'CFF file does not define a top-level version: %s', cffpath);
   elseif numel(matches) > 1
      error('icemodel:internal:readCffVersion:DuplicateVersion', ...
         'CFF file defines more than one top-level version: %s', cffpath);
   end

   % Parse quoted scalars before unquoted comments so a hash inside matching
   % quotes remains part of the version value.
   separator = strfind(matches{1}, ':');
   raw = strtrim(matches{1}(separator(1) + 1:end));
   quoted = false;
   if ~isempty(raw) && raw(1) == ''''
      token = regexp(raw, ...
         '^''((?:[^'']|'''')*)''[ \t]*(?:#.*)?$', 'tokens', 'once');
      if isempty(token)
         error('icemodel:internal:readCffVersion:InvalidVersion', ...
            'CFF version has malformed single quotes: %s', cffpath);
      end
      value = strrep(token{1}, '''''', '''');
      quoted = true;
   elseif ~isempty(raw) && raw(1) == '"'
      token = regexp(raw, '^"([^"]*)"[ \t]*(?:#.*)?$', ...
         'tokens', 'once');
      if isempty(token)
         error('icemodel:internal:readCffVersion:InvalidVersion', ...
            'CFF version has malformed double quotes: %s', cffpath);
      end
      value = token{1};
      quoted = true;
   else
      % YAML starts an unquoted comment at a leading hash or at a hash
      % separated from the scalar by whitespace.
      if ~isempty(raw) && raw(1) == '#'
         value = '';
      else
         comment_start = regexp(raw, '[ \t]+#', 'start', 'once');
         if isempty(comment_start)
            value = raw;
         else
            value = strtrim(raw(1:comment_start - 1));
         end
      end
   end

   % Empty and unquoted YAML null scalars cannot identify an IceModel release;
   % quoted text with the same spelling remains a valid literal version.
   if isempty(value) || (~quoted && ...
         (strcmp(value, '~') || strcmpi(value, 'null')))
      error('icemodel:internal:readCffVersion:InvalidVersion', ...
         'CFF version must not be empty: %s', cffpath);
   end
   version = value;
end
