clean
%------------------------------------------------------------------------------
%   run dependent information
%------------------------------------------------------------------------------
savedata    =   false;
sitenames   =  {'behar'};  
forcingdata =  {'kanm'};   % {'mar','kanl'}
userdata    =  {'none'};   % {'none','merra','racmo','mar','modis'};
uservars    =  {'albedo'};
meltmodels  =  {'icemodel'};
startyear   =  2016;
endyear     =  2016;
simyears    =  startyear:endyear;

ensemble    =  ensembleList(  forcingdata,userdata,uservars,meltmodels, ...
                              simyears,sitenames);

% activate the right version (need to add these functions to the repo)
workoff skinmodel
workon icemodel

%------------------------------------------------------------------------------
% set paths
%------------------------------------------------------------------------------
opts.path.input    = [pwd '/input/'];
opts.path.output   = [pwd '/output/'];
opts.path.metdata  = [opts.path.input 'met/'];
opts.path.initdata = [opts.path.input 'init/'];
opts.path.userdata = [opts.path.input 'userdata/'];

% valid userdata:
% {'racmo','mar','merra','kanm','kanl','modis','none'}

% valid sitenames:
% {'upperBasin','behar','slv1','slv2','ak4'}

% valid forcingdata:
% {'MAR','KANM','KANL'}

%------------------------------------------------------------------------------
% run all combos
%------------------------------------------------------------------------------
for n = 1:ensemble.numcombos
   
   sitename    = char(ensemble.allcombos.sitenames(n));
   simyear     = char(ensemble.allcombos.simyears(n));
   meltmodels  = char(ensemble.allcombos.meltmodels(n));
   forcingdata = char(ensemble.allcombos.forcingdata(n));
   userdata    = char(ensemble.allcombos.userdata(n));
   uservars    = cellstr(ensemble.allcombos.uservars(n));
   
   % this skips cases like 'icemodel_upperbasin_2016_kanl_swap_KANL_albedo'
   % this does not skip the cases where the sitename and forcing are the
   % same, which are the cases where kanm/kanl are used to force basin runs
   if strcmpi(userdata,forcingdata)
      continue;
   end
   
   % set the model options
   opts        = a_opts(   opts,sitename,simyear,meltmodels,forcingdata,...
                           userdata,uservars,str2double(simyear),...
                           str2double(simyear));

   % display the model configuration
   disp([meltmodels ', ' sitename ', '  forcingdata ', ' userdata ', ' simyear])
   
   
%------------------------------------------------------------------------------
%     RUN THE MODEL
%------------------------------------------------------------------------------
   
   if opts.skinmodel == true
      tic; [ice1,ice2,met,opts] = skinmodel(opts); toc
   else
      tic; [ice1,ice2,met,opts] = icemodel(opts); toc
   end

   % post process
   [ice1,ice2,met] = POSTPROC(ice1,ice2,met,opts);

%------------------------------------------------------------------------------
%     save the data
%------------------------------------------------------------------------------
   if savedata == true
      if ~exist(opts.path.output,'dir'); mkdir(opts.path.output); end 
      save([opts.path.output opts.fsave],'ice1','ice2','opts');
   end

end

%------------------------------------------------------------------------------
% quick eval
%------------------------------------------------------------------------------

[Runoff,Discharge,Catchment]  = prepRunoff(opts,ice1);
AblationHourly                = prepAblation(opts,ice1,'hourly');
AblationDaily                 = prepAblation(opts,ice1,'daily');

if strcmpi(userdata,'none')
   ltext = {'ADCP','RACMO','MAR','MERRA',['ICE (' upper(forcingdata) ')']};
else
   ltext = {'ADCP','RACMO','MAR','MERRA',['ICE (' upper(userdata) ')']};
end

% FIGURE 1 = RUNOFF (options: 'raw','mean','members','sensitivity','surf')
if opts.skinmodel == true
	h1 = plotRunoff(Runoff,Discharge,Catchment,  'plotsurf',true,          ...
                                                'legendtext',ltext,     ...
                                                'sitename',sitename,    ...
                                                'refstart',true);
else
   h1 = plotRunoff(Runoff,Discharge,Catchment,'plotensemble',false,  ...
                                            'legendtext',ltext,      ...
                                            'sitename',sitename);
end

% PLOT ABLATION
t1 = datetime(str2double(simyear),6,1,0,0,0,'TimeZone','UTC');
plotPromice(AblationHourly,'refstart',t1);
plotPromice(AblationDaily,'refstart',t1);

% plot enbal - need option to plot met station data when forcing is rcm 
plotEnbal(ice1,met);
