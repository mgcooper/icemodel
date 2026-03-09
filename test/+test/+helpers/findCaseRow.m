function idx = findCaseRow(T, case_id)
%FINDCASEROW Find first row in a table matching CASE_ID.
   arguments
      T table
      case_id
   end

   idx = [];
   if isempty(T) || ~ismember('case_id', T.Properties.VariableNames)
      return
   end
   idx = find(string(T.case_id) == string(case_id), 1, 'first');
end
