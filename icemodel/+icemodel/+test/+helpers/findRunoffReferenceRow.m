function idx = findRunoffReferenceRow(T, c)
   %FINDRUNOFFREFERENCEROW Resolve runoff reference row for one formal case.
   %
   %  idx = icemodel.test.helpers.findRunoffReferenceRow(T, c)
   arguments
      T table
      c
   end

   % Accept either a single table row or an equivalent scalar struct.
   if istable(c)
      assert(height(c) == 1, ...
         'findRunoffReferenceRow expects a single table row')
      c = table2struct(c);
   end

   % Prefer an explicit runoff_site field when the case already resolved it.
   if isfield(c, 'runoff_site') && ~isempty(c.runoff_site)
      sitename = string(c.runoff_site);
   else
      sitename = icemodel.test.helpers.getRunoffSite(c.sitename);
   end

   % Match site, forcing family, and year first, then fall back to site/year.
   idx = find(T.sitename == sitename ...
      & T.forcings == string(c.forcings) ...
      & T.simyear == c.simyear, 1, 'first');

   if isempty(idx)
      idx = find(T.sitename == sitename ...
         & T.simyear == c.simyear, 1, 'first');
   end
end
