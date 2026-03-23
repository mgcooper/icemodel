function idx = findCaseRow(T, case_id)
   %FINDCASEROW Find the first baseline/report row matching CASE_ID.
   %
   %  idx = icemodel.test.helpers.findCaseRow(T, case_id)
   arguments
      T table
      case_id
   end

   % Default to empty so callers can treat "not found" as an empty index.
   idx = [];
   % Treat missing tables or schemas as a clean "not found".
   if isempty(T) || ~ismember('case_id', T.Properties.VariableNames)
      return
   end
   idx = find(string(T.case_id) == string(case_id), 1, 'first');
end
