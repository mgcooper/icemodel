function cases = getPerfCaseMatrix(kwargs)
%GETPERFCASEMATRIX Return the canonical formal performance case matrix.
%
%  cases = icemodel.test.helpers.getPerfCaseMatrix()
%  cases = icemodel.test.helpers.getPerfCaseMatrix(tier="full")
%  cases = icemodel.test.helpers.getPerfCaseMatrix(smbmodel="skinmodel")
%  cases = icemodel.test.helpers.getPerfCaseMatrix(solver=[1 3], simyear=2017)
%  cases = icemodel.test.helpers.getPerfCaseMatrix(smoke_sites="kanm", ...
%     full_sites=["kanm"; "kanl"])
%
% Inputs:
%  tier        - smoke | full | all formal perf subset
%  smbmodel    - all | icemodel | skinmodel
%  solver      - optional subset of [1 2 3]
%  simyear     - single benchmark year used by the perf suite
%  smoke_sites - sites included in the smoke matrix
%  full_sites  - sites included in the full matrix
%
% Output:
%  cases - table with columns:
%    case_id, tier, family, smbmodel, sitename, forcings, userdata, ...
%    uservars, simyear, solver
%
%  `case_id` identifies the underlying model run, not the suite tier. This
%  lets smoke and full compare against the same perf baseline row when they
%  run the same physical case.

   arguments
      kwargs.tier (1, :) string {mustBeMember(kwargs.tier, ...
         ["smoke", "full", "all"])} = "smoke"
      kwargs.smbmodel (1, :) string {mustBeMember(kwargs.smbmodel, ...
         ["all", "icemodel", "skinmodel"])} = "all"
      kwargs.solver {icemodel.validators.mustBeSolverFilter(kwargs.solver)} = []
      kwargs.simyear (1, 1) double {mustBeInteger, mustBePositive} = 2016
      kwargs.smoke_sites string = "kanm"
      kwargs.full_sites string = ["kanm"; "kanl"]
   end
   [tier, smbmodel, solver, simyear, smoke_sites, full_sites] = deal( ...
      kwargs.tier, kwargs.smbmodel, kwargs.solver, kwargs.simyear, ...
      string(kwargs.smoke_sites(:)), string(kwargs.full_sites(:)));

   % Build the fixed smoke/full matrices once from the explicit inputs above.
   smoke = makeCases("smoke", smoke_sites, simyear);
   full = makeCases("full", full_sites, simyear);

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
      cases = cases(ismember(cases.smbmodel, smbmodel), :);
   end

   if ~isempty(solver)
      cases = cases(ismember(cases.solver, solver), :);
   end
end

function cases = makeCases(tier_name, sites, simyear)
   models = icemodel.test.helpers.formalSmbmodels();
   rows = struct([]);
   k = 0;

   for isite = 1:numel(sites)
      sitename = sites(isite);
      for imodel = 1:numel(models)
         smbmodel = models(imodel);
         solver_cases = formalSolversForModel(smbmodel);
         for isolver = 1:numel(solver_cases)
            solver_id = solver_cases(isolver);
            k = k + 1;
            rows(k).case_id = icemodel.test.helpers.makeFormalCaseId( ...
               smbmodel, sitename, simyear, solver_id);
            rows(k).tier = string(tier_name);
            rows(k).family = "self";
            rows(k).smbmodel = smbmodel;
            rows(k).sitename = sitename;
            rows(k).forcings = sitename;
            rows(k).userdata = "";
            rows(k).uservars = "";
            rows(k).simyear = simyear;
            rows(k).solver = solver_id;
         end
      end
   end

   cases = struct2table(rows);
end

function solver_cases = formalSolversForModel(smbmodel)
   if smbmodel == "icemodel"
      solver_cases = 1:3;
   else
      solver_cases = 1;
   end
end
