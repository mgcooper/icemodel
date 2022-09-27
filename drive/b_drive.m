clean
%------------------------------------------------------------------------------
%   run dependent information
%------------------------------------------------------------------------------
ID          =  'run1'; % run1 = w/chi, run2 = w/o chi (chi=0)
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

% activate the right version
workoff skinmodel
workon icemodel

% note: ak4 userData merra and racmo is available from 2012:2016, whereas
% mar, kanl, and modis are available from 2009:2016

%------------------------------------------------------------------------------
% set paths
%------------------------------------------------------------------------------
opts.path.input    = setpath('GREENLAND/icemodel/input/');
opts.path.output   = ['/Users/coop558/mydata/icemodel/output/v10/' ID '/'];
opts.path.metdata  = [opts.path.input 'met/'];
opts.path.initdata = [opts.path.input 'init/'];
opts.path.userdata = [opts.path.input 'userdata/'];

% valid userdata:
% {'racmo','mar','merra','kanm','kanl','modis','none'}

% valid sitenames:
% {'upperBasin','behar','slv1','slv2','ak4'}

% valid forcingdata:
% {'MAR','KANM','KANL'}

% MAR forcing + MAR, MERRA, RACMO, AWS, MODIS albedo
% MERRA forcing + MAR, MERRA, RACMO, AWS, MODIS albedo
% RACMO forcing + MAR, MERRA, RACMO, AWS, MODIS albedo
% AWS forcing + MAR, MERRA, RACMO, AWS, MODIS albedo

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
   
   % this skips cases like 'icemodel_upperBasin_2016_kanl_swap_KANL_albedo'
   % this does not skip the cases where the sitename and forcing are the
   % same, which are the cases where kanm/kanl are used to force basin runs
   if strcmpi(userdata,forcingdata)
      continue;
   end
   
   % set the model options (including the fsave prefix)
   opts        = a_opts(   opts,sitename,simyear,meltmodels,forcingdata,...
                           userdata,uservars,str2double(simyear),...
                           str2double(simyear));

   % run the model 

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
   
% % if using diags:   
%    tic; [ice1,ice2,met,opts,diags] = icemodel(opts); toc   
%    [ice1,ice2,met,diags] = POSTPROC(ice1,ice2,met,opts,diags);


%------------------------------------------------------------------------------
%     save the data
%------------------------------------------------------------------------------
   if savedata == true
      if ~exist(opts.path.output,'dir'); mkdir(opts.path.output); end 
      save([opts.path.output opts.fsave],'enbal','ice1','ice2','opts');
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
   if ismember(sitename,{'slv1','slv2'})
      h1 = plotRunoff(Runoff,Discharge,Catchment,'plotensemble',false,  ...
                                                'legendtext',ltext,     ...
                                                'sitename',sitename,    ...
                                                'refstart',true);
      h1 = plotRunoff(Runoff,Discharge,Catchment,'plotensemble',false,  ...
                                                'legendtext',ltext,     ...
                                                'sitename',sitename,    ...
                                                'refstart',false);
   else
      h1 = plotRunoff(Runoff,Discharge,Catchment,'plotensemble',false,  ...
                                               'legendtext',ltext,      ...
                                               'sitename',sitename);
   end
end

% PLOT ABLATION
t1 = datetime(str2double(simyear),6,1,0,0,0,'TimeZone','UTC');
h2 = plotPromice(AblationHourly,'refstart',t1);
h2 = plotPromice(AblationDaily,'refstart',t1);

t1 = datetime(str2double(simyear),7,1,0,0,0,'TimeZone','UTC');
h2 = plotPromice(AblationHourly,'refstart',t1);
h2 = plotPromice(AblationDaily,'refstart',t1);

% plot enbal - need option to plot met station data when forcing is rcm 
plotEnbal(ice1,met);

% plotTice(opts,ice2);
% 
% % if mar forcing, then met.tsfc will be daily
% 
% t1 =  datetime(year(ice1.Time(1)),1,1,0,0,0,'TimeZone','UTC');
% t2 =  datetime(year(ice1.Time(1)),12,31,23,0,0,'TimeZone','UTC');
% 
% plotice2(ice2,'Tice','t1',t1,'t2',t2);


% this is the percent of timesteps that Tsfc did not converge
% 100*(1-(sum(diags.Tflag)/opts.maxiter))
% 
% figure; plot(enbal.tsfc,enbal.TN,'o'); addOnetoOne; xylabel('Tsfc','TN');
% set(gca,'XLim',[-10 1],'YLim',[-10 1]);
% 
% % want to know Qc for timesteps when Tsfc does not converge
% figure; plot(diags.Time,diags.Tflag)
% 
% mean(diags.Qc(diags.Tflag==true))
% mean(diags.Qc(diags.Tflag==false))
% 
% max(diags.Qc(diags.Tflag==true))
% max(diags.Qc(diags.Tflag==false))
% 
% min(diags.Qc(diags.Tflag==true))
% min(diags.Qc(diags.Tflag==false))
% 
% histogram(diags.Qc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pause; close all;

% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 
% % not 100% sure but I think this converts the W/m2 enbalance error to m.w.e
% load( 'PHYSCONS', 'cp_ice','cp_liq','Lf');
% cp_sno   = interp1(met.Time,ice2.cp_sno(1,:),ice1.Time);
% werr     = cp_sno.*enbal.balance.*opts.dz_thermal./Lf;
% % figure; plot(enbal.Time,cumsum(werr)); ylabel('cumulative error (m w.e.)')


ice1.runoff(end)
try
   printf(h1.data.("ICE (MODIS)")(end),0)
   printf((h1.data.("ICE (MODIS)")(end)-h1.data.ADCP(end))/h1.data.ADCP(end),2)
end

try
   printf(h1.data.("ICE (KANM)")(end),0)
   printf((h1.data.("ICE (KANM)")(end)-h1.data.ADCP(end))/h1.data.ADCP(end),2)
end

try
   printf(h1.data.("ICE (NONE)")(end),0)
   printf((h1.data.("ICE (NONE)")(end)-h1.data.ADCP(end))/h1.data.ADCP(end),2)
end

try
   printf(h1.data.("ICE (KANL)")(end),0)
   printf((h1.data.("ICE (KANL)")(end)-h1.data.ADCP(end))/h1.data.ADCP(end),2)
end

% % note: promice ice ablation: 
% % kanL 2015 = 
% % kanM 2015 = 0 (Fausto)
% % kanM 2015 = 0 (Me)
% % kanM 2016 = 2.3 (Fausto)
% % kanM 2016 = 2.1 (Me, m.w.e) 2.3 (Me, m.i.e)


% % figure; 
% % plot(ice1.Time,ice1.MFerr); hold on;
% % plot(ice1.Time,ice1.LCerr); 
% % legend('MF error','LC error');
% % ylabel('error (^oC)');



% % this doesn't work quite as hoped b/c of the situation where kanl/kanm
% % are swapped for mar in the forcing data step of metiinit
% %    % if the met file does not exist, continue
%    fmet  = [opts.path.metData opts.metfname];
%    if exist(fmet,'file') ~= 2 ; continue; end
% % 
% %    % activating this option prevents re-running a site so use with care
% %    % or the data file already exists
% %    fsave = [opts.path.output opts.fsave '.mat'];
% %    if exist(fsave,'file') == 2 ; continue; end
% 
% % check if the forcingUserData file exists
%    userData    =   opts.forcingUserData;
%    userPath    =   opts.path.userData;
%    userFile    =   [userPath userData '_' siteName '_' simYear '.mat']; 
%    if exist(userFile,'file') ~= 2 ; continue; end
