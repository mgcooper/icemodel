function mustBeStationToken(tokens)
   %MUSTBESTATIONTOKEN Require canonical lowercase alphanumeric station IDs.
   %
   %  Public reconstruction and report boundaries call this validator before
   %  station tokens enter globs or output paths. A caller must handle the
   %  sentinel values "auto" and "all" before it calls this validator.

   tokens = string(tokens);
   valid = arrayfun(@(token) ~isempty(regexp(char(token), ...
      '^[a-z][a-z0-9]*$', 'once')), tokens);
   if ~all(valid, 'all')
      error('icemodel:reconstruct:mustBeStationToken:invalidToken', ...
         'station tokens must match ^[a-z][a-z0-9]*$');
   end
end
