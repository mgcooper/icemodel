function mustBeStationToken(tokens)
   %MUSTBESTATIONTOKEN Require canonical lowercase alphanumeric station IDs.
   %
   %  Used at public reconstruction/report boundaries before station tokens
   %  enter globs or output paths. Sentinel values such as "auto" and "all"
   %  must be handled by callers before invoking this validator.

   tokens = string(tokens);
   valid = arrayfun(@(token) ~isempty(regexp(char(token), ...
      '^[a-z][a-z0-9]*$', 'once')), tokens);
   if ~all(valid, 'all')
      error('icemodel:reconstruct:mustBeStationToken:invalidToken', ...
         'station tokens must match ^[a-z][a-z0-9]*$');
   end
end
