clearvars
close all
clc
%#ok<*ASGLU>

% This script demonstrates how to run icemodel (or skinmodel) at a point. Note
% that the examples below are set up to run skinmodel (the surface energy
% balance model) because it runs fast. Replace "skinmodel" with "icemodel" to
% run icemodel, which includes the subsurface energy balance.
%
% IceModel requires a met (forcing) file with the following naming convention:
% ['met_' sitename '_' forcings '_' yearname '_' timestep '.mat']
%
% The example file included with this repo: 'met_KANM_KANM_2016_1hr.mat'
%
% This met file contains a timetable with forcings from the KANM weather
% station for year 2016 on a 1 hr timestep. Note the sitename and the
% forcings are the same in this example.
%
% To run icemodel (or skinmodel) for your site, create a metfile with the same
% format as the example one, save it in the input/met folder, and add the
% sitename to the namelists in the icemodel.run.point arguments block.

% Source the project-level configuration.
config = icemodel.config(casename="demo");

%% Example 1: KAN-M weather station in 2016

% Set the run configuration
saveflag = false;
sitename = 'kanm';
forcings = 'kanm';
smbmodel = 'skinmodel';
simyears = 2016:2016;
backupflag = false;

% Run the model
[ice1, ice2, met, opts] = icemodel.run.point(...
   saveflag=saveflag, ...
   sitename=sitename, ...
   forcings=forcings, ...
   smbmodel=smbmodel, ...
   simyears=simyears, ...
   backupflag=backupflag);

% Demo plot
icemodel.plot.enbal(ice1, met)

%% Example 2 - Replace KAN-M albedo with MODIS albedo

% Set the run configuration. Use the "userdata" option to swap the albedo.
saveflag = false;
sitename = 'kanm';
forcings = 'kanm';
smbmodel = 'skinmodel';
userdata = 'modis';
uservars = 'albedo';
simyears = 2016:2016;
backupflag = false;

% Run the model
[ice1, ice2, met, opts] = icemodel.run.point(...
   saveflag=saveflag, ...
   sitename=sitename, ...
   forcings=forcings, ...
   smbmodel=smbmodel, ...
   simyears=simyears, ...
   userdata=userdata, ...
   uservars=uservars, ...
   backupflag=backupflag);

% Demo plot
icemodel.plot.enbal(ice1, met)

%% Example 3 - Replace KAN-M air temperature with MERRA-2 air temperature

% Set the run configuration. Use the "userdata" option to swap the air temp.
saveflag = false;
sitename = 'kanm';
forcings = 'kanm';
smbmodel = 'skinmodel';
userdata = 'merra';
uservars = 'tair';
simyears = 2016:2016;
backupflag = false;

% Run the model
[ice1, ice2, met, opts] = icemodel.run.point(...
   saveflag=saveflag, ...
   sitename=sitename, ...
   forcings=forcings, ...
   smbmodel=smbmodel, ...
   simyears=simyears, ...
   userdata=userdata, ...
   uservars=uservars, ...
   backupflag=backupflag);

% Demo plot
icemodel.plot.enbal(ice1, met)

