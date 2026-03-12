function results = run_regression_suite(kwargs)
   %RUN_REGRESSION_SUITE Run formal icemodel numerical regression test suite.
   %
   %  results = run_regression_suite()
   %  results = run_regression_suite(tier="smoke")
   %  results = run_regression_suite(tier="full")
   %  results = run_regression_suite(tier="smoke", smbmodel="skinmodel")
   %  results = run_regression_suite(tier="smoke", smbmodel="icemodel", solver=2)
   %  results = run_regression_suite(tier="smoke", smbmodel="icemodel", solver=[2 3])
   %  results = run_regression_suite(tier="smoke", smbmodel="icemodel", solver=[1 3])
   %  results = run_regression_suite(tier="full", baseline="v1.1")
   %
   % Use this for normal regression comparisons against an existing rolling
   % or release baseline. This function does not update baselines; it only
   % runs the formal cases, compares core scalar outputs to the requested
   % baseline, and writes one artifact under test/artifacts/<run_name>/.
   % The optional solver filter accepts any subset of [1 2 3].
   %
   % CLI entrypoint:
   %  matlab -batch "run('/ABS/PATH/icemodel/test/run_regression_suite.m')"

   arguments (Input)
      kwargs.tier (1, :) string {mustBeMember(kwargs.tier, ...
         ["smoke", "full", "all"])} = "smoke"
      kwargs.smbmodel (1, :) string {mustBeMember(kwargs.smbmodel, ...
         ["all", "icemodel", "skinmodel"])} = "all"
      kwargs.solver {mustBeValidSolverFilter(kwargs.solver)} = []
      kwargs.baseline string = string.empty()
      kwargs.run_name string = string.empty()
   end
   [tier, smbmodel, solver, baseline, run_name] = deal(kwargs.tier, ...
      kwargs.smbmodel, kwargs.solver, kwargs.baseline, kwargs.run_name);

   % Provide the requested regression selection to the unittest class.
   setenv('ICEMODEL_TEST_TIER', char(tier));
   setenv('ICEMODEL_TEST_SMBMODEL_FILTER', char(smbmodel));
   setenv('ICEMODEL_TEST_SOLVER_FILTER', char(join(string(solver), ',')));
   setenv('ICEMODEL_REGRESSION_BASELINE', char(baseline));
   setenv('ICEMODEL_TEST_RUN_NAME', char(run_name));

   % Ensure this test folder is on path if called from CLI.
   thisdir = fileparts(mfilename('fullpath'));
   addpath(thisdir);

   % Run the unittest class with standard MATLAB text output.
   suite = testsuite(fullfile(thisdir, 'IcemodelRegressionTest.m'));
   runner = matlab.unittest.TestRunner.withTextOutput;
   results = runner.run(suite);
end

function mustBeValidSolverFilter(x)
   if isempty(x)
      return
   end
   mustBeNumeric(x)
   mustBeInteger(x)
   if any(~ismember(x, [1 2 3]))
      error('solver must be empty or a subset of [1 2 3]')
   end
end
