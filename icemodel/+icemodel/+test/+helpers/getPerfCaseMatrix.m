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
   %  smoke_sites - advanced override for the smoke-tier site list
   %  full_sites  - advanced override for the full-tier site list
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
      kwargs.tier (1, :) string ...
         {icemodel.validators.mustBeTestTierName(kwargs.tier)} = "smoke"
      kwargs.smbmodel (1, :) string ...
         {icemodel.validators.mustBeTestSmbmodelSelector(kwargs.smbmodel)} = "all"
      kwargs.solver {icemodel.validators.mustBeSolverFilter(kwargs.solver)} = []
      kwargs.simyear (1, 1) double {mustBeInteger, mustBePositive} = 2016
      kwargs.smoke_sites string = "kanm"
      kwargs.full_sites string = ["kanm"; "kanl"]
   end

   % Deal out arguments.
   [tier, smbmodel, solver, simyear, smoke_sites, full_sites] = deal( ...
      kwargs.tier, kwargs.smbmodel, kwargs.solver, kwargs.simyear, ...
      reshape(kwargs.smoke_sites, [], 1), reshape(kwargs.full_sites, [], 1));

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

   % Apply the explicit model and solver filters after constructing the
   % stable smoke/full matrices.
   if ~any(strcmpi(smbmodel, "all"))
      cases = cases(ismember(cases.smbmodel, smbmodel), :);
   end

   if ~isempty(solver)
      % The solver filter only narrows the icemodel rows. Keep single-solver
      % models such as skinmodel in the matrix when smbmodel="all" so the
      % shared perf suite still covers the full model family.
      is_icemodel = cases.smbmodel == "icemodel";
      keep = ~is_icemodel | ismember(cases.solver, solver);
      cases = cases(keep, :);
   end
end

function cases = makeCases(tier_name, sites, simyear)
   %MAKECASES Expand one tier/site selection into the formal perf case rows.
   models = icemodel.namelists.smbmodel("test");
   rows = struct([]);
   k = 0;

   % Each site/model pair contributes one row per supported solver.
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
   %FORMALSOLVERSFORMODEL Return the solver ids covered by the perf suite.
   if smbmodel == "icemodel"
      solver_cases = 1:3;
   else
      solver_cases = 1;
   end
end
