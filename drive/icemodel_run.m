%clean

% this demonstrates how to set up a one-year run at one site

%------------------------------------------------------------------------------
%  set the run-specific model configuration
%------------------------------------------------------------------------------
savedata    = false;

% these settings are for the rio behar catchment:

sitename    = 'behar';        % options: 'kanm', 'behar'
startyear   = 2016;           % options: 2016
endyear     = 2016;           % options: 2016
meltmodel   = 'icemodel';     % options: 'icemodel','skinmodel'
forcingdata = 'kanm';         % options: 'mar','kanm'
userdata    = 'modis';        % options: 'modis','racmo','merra','mar','kanm','none'
uservars    = 'albedo';       % options: 'albedo', or any var in met

%------------------------------------------------------------------------------
%  set the input and output paths (see icemodel_config.m)
%------------------------------------------------------------------------------
icemodel_config

%------------------------------------------------------------------------------
%  build the 'opts' model configuration structure
%------------------------------------------------------------------------------
opts = icemodel_opts(   sitename,meltmodel,forcingdata,userdata,uservars,  ...
                        startyear,endyear);

%------------------------------------------------------------------------------
%  run the model
%------------------------------------------------------------------------------
[ice1,ice2,met,opts] = icemodel(opts);

%------------------------------------------------------------------------------
%  save the data
%------------------------------------------------------------------------------
if savedata == true
   if ~exist(opts.pathoutput,'dir'); mkdir(opts.pathoutput); end
   save([opts.pathoutput opts.fsave],'ice1','ice2','opts');
end

%------------------------------------------------------------------------------
%  prep the output for plotting
%------------------------------------------------------------------------------
[Runoff,Discharge,Catchment]  = prepRunoff(opts,ice1);
AblationHourly                = prepAblation(opts,ice1,'hourly');
AblationDaily                 = prepAblation(opts,ice1,'daily');

t1 = datetime(startyear,7,1,0,0,0,'TimeZone','UTC');
t2 = datetime(startyear,9,1,0,0,0,'TimeZone','UTC');

% plot the runoff
if opts.meltmodel == "skinmodel"
	h1 = plotRunoff(Runoff,Discharge,Catchment,'plotsurf',true,'sitename',  ...
                     sitename,'userdata',userdata,'forcingdata',forcingdata);
else
   h1 = plotRunoff(Runoff,Discharge,Catchment,'sitename',sitename,         ...
         'userdata',userdata,'forcingdata',forcingdata,'t1',t1,'t2',t2);
end

% plot ablation
t1 = datetime(startyear,6,1,0,0,0,'TimeZone','UTC');
plotPromice(AblationDaily,'refstart',t1);
plotPromice(AblationHourly,'refstart',t1);

% plot the energy balance
% plotEnbal(ice1,met);


