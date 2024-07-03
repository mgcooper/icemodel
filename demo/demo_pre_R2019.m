clearvars
close all
clc
%#ok<*ASGLU>

% This script demonstrates how to run icemodel (or skinmodel) at a point, in a
% way that is compatible with MATLAB versions back to at least R2017a.

% Run this script from the top-level icemodel repo folder.

%% Set the project paths

% Run setup() to add the project paths to the matlab path. Note that setup()
% also sets the ICEMODEL_INPUT_PATH and ICEMODEL_OUTPUT_PATH environment
% variables to the demo/input and demo/output folders.
setup()

%% Example 1: KAN-M weather station in 2016

% Set the run configuration.
saveflag = false;
sitename = 'kanm';
forcings = 'kanm';
smbmodel = 'skinmodel';
simyears = 2016:2016;
backupflag = false;

% Keep these options empty for this demonstration.
userdata = [];
uservars = [];
testname = [];

% To run the model in a way that is compatible with matlab versions pre-2019a,
% directly call icemodel.setopts and pass the options to the model function.
% Note that for >=R2019a, the convenience function icemodel.run.point can be
% used to perform each of these steps, with additional type checking via the
% arguments block.

% Set the model options.
opts = icemodel.setopts(smbmodel, sitename, simyears, forcings, ...
   userdata, uservars, testname, saveflag, backupflag);

% Pass the model options to the skinmodel (or icemodel) function
switch smbmodel
   case 'icemodel'
      tic; [ice1, ice2] = icemodel(opts); toc
   case 'skinmodel'
      tic; [ice1, ice2] = skinmodel(opts); toc
end

% Run the post-processing step.
[ice1, ice2, met] = POSTPROC(ice1, ice2, opts, simyears);

% Create the demo plot of the surface energy balance compared with observed
% forcings from the KAN-M weather station.
icemodel.plot.enbal(ice1, met)
