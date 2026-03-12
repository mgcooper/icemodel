function case_id = makeFormalCaseId(smbmodel, sitename, simyear, solver)
%MAKEFORMALCASEID Return canonical formal-suite identifier for one model run.
%
%  case_id = test.helpers.makeFormalCaseId("icemodel", "kanm", 2016, 2)
%
% The identifier is intentionally independent of suite tier (smoke/full) and
% suite type (perf/regression). It represents the underlying physical run so
% the same baseline row can be reused across different formal comparisons.
   arguments
      smbmodel (1, 1) string
      sitename (1, 1) string
      simyear (1, 1) double {mustBeInteger, mustBePositive}
      solver (1, 1) double {mustBeInteger, mustBePositive}
   end

   case_id = smbmodel + "_" + sitename + "_" + string(simyear) ...
      + "_solver" + string(solver);
end
