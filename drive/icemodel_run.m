
% this demonstrates how to set up a single-year run at one site

% set the main configuration options
sitename    = 'KANM';                  % KAN-M weather station
startyear   = 2018;
endyear     = 2018;
meltmodel   = 'icemodel';
forcingdata = 'KANM';
userdata    = 'modis';                 % use modis albedo instead of kan-m
uservars    = 'albedo';

% set the input and output paths (see icemodel_config.m)
setenv('ICEMODELIINPUTPATH','/full/path/to/icemodel/input/');
setenv('ICEMODELOUTPUTPATH','/full/path/to/icemodel/output/');

% set the model options 
opts = icemodel_opts(sitename,meltmodel,forcingdata,userdata,uservars,startyear,endyear);

% run the model
[ice1,ice2,met,opts] = icemodel(opts);

% save the data

% plot the data