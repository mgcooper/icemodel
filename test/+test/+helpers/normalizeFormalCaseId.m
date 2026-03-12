function case_id = normalizeFormalCaseId(case_id)
%NORMALIZEFORMALCASEID Normalize legacy formal-suite identifiers.
%
%  case_id = test.helpers.normalizeFormalCaseId(case_id)
%
% Older saved perf/regression baselines encoded suite prefixes (smoke_/full_/
% reg_) and used `_bcN` for the solver index. Normalize those forms so current
% loaders can still match legacy baseline rows.
   case_id = string(case_id);
   for i = 1:numel(case_id)
      s = case_id(i);
      s = regexprep(s, "^(smoke|full|reg)_", "");
      tok = regexp(s, "^(.*)_bc([0-9]+)$", "tokens", "once");
      if ~isempty(tok)
         s = string(tok{1}) + "_solver" + string(tok{2});
      end
      case_id(i) = s;
   end
end
