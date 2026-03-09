function idx = findRunoffReferenceRow(T, c)
%FINDRUNOFFREFERENCEROW Resolve runoff reference row for one formal case.
   arguments
      T table
      c
   end

   if istable(c)
      assert(height(c) == 1, ...
         'findRunoffReferenceRow expects a single table row')
      c = table2struct(c);
   end

   if isfield(c, 'runoff_site') && ~isempty(c.runoff_site)
      sitename = string(c.runoff_site);
   else
      sitename = test.helpers.getRunoffSite(c.sitename);
   end

   idx = find(T.sitename == sitename ...
      & T.forcings == string(c.forcings) ...
      & T.simyear == c.simyear, 1, 'first');

   if isempty(idx)
      idx = find(T.sitename == sitename ...
         & T.simyear == c.simyear, 1, 'first');
   end
end
