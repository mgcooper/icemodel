function cases = getRegressionCaseMatrix(tier, smbmodel, solver)
%GETREGRESSIONCASEMATRIX Return deterministic regression test case matrix.
%
%  cases = test.helpers.getRegressionCaseMatrix("smoke")
%  cases = test.helpers.getRegressionCaseMatrix("full")
%  cases = test.helpers.getRegressionCaseMatrix("all")
%  cases = test.helpers.getRegressionCaseMatrix("smoke", "skinmodel")
%  cases = test.helpers.getRegressionCaseMatrix("smoke", "icemodel", 2)
%
% The formal regression matrix is intentionally compact and stable:
% self-forced station runs at `kanm` and `kanl` for 2016. Icemodel is crossed
% with `solver = 1:3`, while skinmodel uses its default path with
% `solver = 1`.
   arguments
      tier = "smoke"
      smbmodel = "all"
      solver = []
   end

   tier = string(tier);
   smbmodel = string(smbmodel);

   % Keep the formal regression matrix compact and fixed for stable baselines.
   smoke = makeCases("regression_smoke", "kanm");
   full = makeCases("regression_full", ["kanm"; "kanl"]);

   switch lower(char(tier))
      case 'smoke'
         cases = smoke;
      case 'full'
         cases = full;
      case 'all'
         cases = [smoke; full];
      otherwise
         error('unrecognized regression tier: %s (expected smoke|full|all)', tier)
   end

   if ~any(strcmpi(smbmodel, "all"))
      % Filter after construction so one helper defines the canonical cases.
      cases = cases(ismember(cases.smbmodel, smbmodel), :);
   end

   if ~isempty(solver)
      cases = cases(ismember(cases.solver, solver), :);
   end
end

function cases = makeCases(tier_name, sites)
   year = 2016;
   models = ["icemodel"; "skinmodel"];
   n = numel(sites) * 4;

   case_id = strings(n, 1);
   tier = strings(n, 1);
   family = strings(n, 1);
   smbmodel = strings(n, 1);
   sitename = strings(n, 1);
   forcings = strings(n, 1);
   userdata = strings(n, 1);
   uservars = strings(n, 1);
   simyear = zeros(n, 1);
   solver = zeros(n, 1);
   runoff_site = strings(n, 1);

   k = 0;
   for i = 1:numel(sites)
      for im = 1:numel(models)
         if models(im) == "icemodel"
            bcs = 1:3;
         else
            bcs = 1;
         end
         for j = 1:numel(bcs)
            k = k + 1;
            simyear(k) = year;
            case_id(k) = test.helpers.makeFormalCaseId( ...
               models(im), sites(i), simyear(k), bcs(j));
            tier(k) = tier_name;
            family(k) = "self";
            smbmodel(k) = models(im);
            sitename(k) = sites(i);
            forcings(k) = sites(i);
            userdata(k) = "";
            uservars(k) = "";
            solver(k) = bcs(j);
            runoff_site(k) = test.helpers.getRunoffSite(sites(i));
         end
      end
   end

   cases = table(case_id, tier, family, smbmodel, sitename, forcings, ...
      userdata, uservars, simyear, solver, runoff_site);
end
