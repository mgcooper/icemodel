function cases = getCaseMatrix(tier, smbmodel)
%GETCASEMATRIX Return deterministic performance test case matrix.
%
%  cases = test.helpers.getCaseMatrix("smoke")
%  cases = test.helpers.getCaseMatrix("full")
%  cases = test.helpers.getCaseMatrix("all")
%  cases = test.helpers.getCaseMatrix("smoke", "skinmodel")
%
% Output:
%  cases - table with columns:
%    case_id, tier, family, smbmodel, sitename, forcings, userdata, uservars,
%    simyear, solver
   arguments
      tier = "smoke"
      smbmodel = "all"
   end

   tier = string(tier);
   smbmodel = string(smbmodel);

   % Build the small fixed smoke matrix and the broader full matrix once.
   smoke = makeSmokeCases();
   full = makeFullCases();

   switch lower(char(tier))
      case 'smoke'
         cases = smoke;
      case 'full'
         cases = full;
      case 'all'
         cases = [smoke; full];
      otherwise
         error('unrecognized tier: %s (expected smoke|full|all)', tier)
   end

   if ~any(strcmpi(smbmodel, "all"))
      % Filter after construction so one helper defines the canonical cases.
      cases = cases(ismember(cases.smbmodel, smbmodel), :);
   end
end

function T = makeSmokeCases()
   bc = [(1:3)'; 1];
   smbmodel = [repmat("icemodel", 3, 1); "skinmodel"];
   n = numel(bc);
   case_id = strings(n, 1);
   for i = 1:n
      case_id(i) = "smoke_" + smbmodel(i) + "_kanm_2016_bc" + int2str(bc(i));
   end

   T = table( ...
      case_id, ...
      repmat("smoke", n, 1), ...
      repmat("self", n, 1), ...
      smbmodel, ...
      repmat("kanm", n, 1), ...
      repmat("kanm", n, 1), ...
      repmat("", n, 1), ...
      repmat("", n, 1), ...
      repmat(2016, n, 1), ...
      bc, ...
      'VariableNames', {'case_id', 'tier', 'family', 'smbmodel', 'sitename', ...
      'forcings', 'userdata', 'uservars', 'simyear', 'solver'});
end

function T = makeFullCases()
   sites = ["kanm"; "kanl"];
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

   k = 0;
   for is = 1:numel(sites)
      s = sites(is);
      for im = 1:numel(models)
         if models(im) == "icemodel"
            bcs = 1:3;
         else
            bcs = 1;
         end
         for ibc = 1:numel(bcs)
            b = bcs(ibc);
            k = k + 1;
            tier(k) = "full";
            family(k) = "self";
            smbmodel(k) = models(im);
            sitename(k) = s;
            forcings(k) = s;
            userdata(k) = "";
            uservars(k) = "";
            simyear(k) = 2016;
            solver(k) = b;
            case_id(k) = "full_" + models(im) + "_" + s + "_2016_bc" + int2str(b);
         end
      end
   end

   T = table(case_id, tier, family, smbmodel, sitename, forcings, userdata, ...
      uservars, simyear, solver);
end
