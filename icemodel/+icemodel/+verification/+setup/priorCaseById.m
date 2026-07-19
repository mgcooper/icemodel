function entry = priorCaseById(cases, case_id)
   %PRIORCASEBYID Return one prior manifest case by canonical case id.
   %
   %  entry = icemodel.verification.setup.priorCaseById(cases, case_id)

   % Missing or legacy case arrays have no reusable identity contract.
   entry = struct([]);
   if isempty(cases) || ~isfield(cases, 'case_id')
      return
   end

   % Preserve the manifest entry verbatim; callers decide which fields refresh.
   hit = find(string({cases.case_id}) == string(case_id), 1);
   if ~isempty(hit)
      entry = cases(hit);
   end
end
