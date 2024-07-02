classdef IcemodelPerfTest < matlab.perftest.TestCase

%    properties
%       opts
%       saveflag = false;
%       sitename = 'kanm';
%       forcings = 'kanm';
%       userdata = [];
%       uservars = 'albedo';
%       smbmodel = 'skinmodel';
%       simyears = 2016;
%    end
% 
%    properties (TestParameter)
%       solver = {1, 2, 3}
%    end
% 
%    methods(TestMethodSetup)
%       function setOpts(testCase)
% 
%          testCase.opts = icemodel.setopts( ...
%             testCase.smbmodel, ...
%             testCase.sitename, ...
%             testCase.simyears, ...
%             testCase.forcings, ...
%             testCase.userdata, ...
%             testCase.uservars, ...
%             testCase.saveflag);
%       end
%    end
% 
%    methods(Test)
%       function testSolver(testCase, solver)
% 
%          testCase.opts.solver = solver;
% 
%          % standard approach (for "slow" code)
%          testCase.startMeasuring
%          if testCase.smbmodel == "icemodel"
%             icemodel(testCase.opts);
%          elseif testCase.smbmodel == "skinmodel"
%             skinmodel(testCase.opts);
%          end
%          testCase.stopMeasuring
%       end
%    end
end
