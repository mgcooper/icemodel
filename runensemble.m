clean

savedata = true;

runcombo = 5;

% note: ak4 userdata merra and racmo is available from 2012:2016, whereas
% mar, kanl, and modis are available from 2009:2016

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

%% Set which predefined
switch runcombo
   case 1
      % DONE upperbasin, 2016, mar, kanl, icemodel
      sitename = {'upperbasin'};
      forcings = {'mar', 'kanl'};
      userdata = {'none','merra','racmo','mar','modis', 'kanl'};
      uservars = {'albedo'};
      simmodel = {'icemodel', 'skinmodel'};
      simyears = 2016:2016;

   case 2
      % DONE behar, 2016, mar, kanm
      sitename = {'behar'};
      forcings = {'mar', 'kanm'};
      userdata = {'none','merra','racmo','mar','modis', 'kanm'};
      uservars = {'albedo'};
      simmodel = {'icemodel', 'skinmodel'};
      simyears = 2016:2016;

   case 3
      % DONE behar, slv1, slv2, 2015, mar, kanm, icemodel
      sitename = {'behar', 'slv1', 'slv2'};
      forcings = {'mar', 'kanm'};
      userdata = {'none','merra','racmo','mar','modis', 'kanm'};
      uservars = {'albedo'};
      simmodel = {'icemodel', 'skinmodel'};
      simyears = 2015:2015;

   case 4
      % DONE ak4, 2009-2011, mar, kanl, icemodel
      sitename = {'ak4'};
      forcings = {'mar', 'kanl'};
      userdata = {'none','mar','modis','kanl'};
      uservars = {'albedo'};
      simmodel = {'icemodel', 'skinmodel'};
      simyears = 2009:2011;

   case 5
      % DONE ak4, 2012-2016, mar, kanl, icemodel, merra, racmo
      sitename = {'ak4'};
      forcings = {'mar', 'kanl'};
      userdata = {'none','merra','racmo','mar','modis','kanl'};
      uservars = {'albedo'};
      simmodel = {'icemodel', 'skinmodel'};
      simyears = 2012:2016;

   otherwise
      % Custom configuration
      sitename = {'ak4'};
      forcings = {'mar', 'kanl'};
      userdata = {'none','merra','racmo','mar','modis','kanl'};
      uservars = {'albedo'};
      simmodel = {'icemodel', 'skinmodel'};
      simyears = 2012:2016;

end

%% Build the ensemble list
ensemble = ensembleList( ...
   forcings, userdata, uservars, simmodel, simyears, sitename);

%% Run all combos
for n = 1:ensemble.numcombos

   sitename = char(ensemble.allcombos.sitename(n));
   simyears = char(ensemble.allcombos.simyears(n));
   simmodel = char(ensemble.allcombos.simmodel(n));
   forcings = char(ensemble.allcombos.forcings(n));
   userdata = char(ensemble.allcombos.userdata(n));
   uservars = char(ensemble.allcombos.uservars(n));

   % this skips cases like 'icemodel_upperBasin_2016_kanl_swap_KANL_albedo'
   % this does not skip the cases where the sitename and forcing are the
   % same, which are the cases where kanm/kanl are used to force basin runs
   if strcmpi(userdata, forcings)
      continue;
   end

   % set the model options
   opts = icemodel.setopts(simmodel, sitename, str2double(simyears), ...
      forcings, userdata, uservars, savedata);

   % display the run info
   disp([simmodel ', ' sitename ', '  forcings ', ' userdata ', ' simyears])

   % RUN THE MODEL
   switch simmodel
      case 'icemodel'
         tic; [ice1, ice2] = icemodel(opts); toc
      case 'skinmodel'
         tic; [ice1, ice2] = skinmodel(opts); toc
   end
end

% post process
met = icemodel.loadmet(opts, numel(str2double(simyears)));
[ice1, ice2] = POSTPROC(ice1, ice2, opts, ...
   met.swd, met.lwd, met.albedo, met.Time);
met = icemodel.processmet(met);
plotice2(ice2, 'Tice')


% quick eval
[Runoff, Discharge, Catchment] = prepRunoff(opts, ice1);
AblationHourly = prepAblation(opts, ice1, 'hourly');
AblationDaily = prepAblation(opts, ice1, 'daily');

if strcmpi(userdata, 'none')
   ltext = {'ADCP','RACMO','MAR','MERRA',['ICE (' upper(forcings) ')']};
else
   ltext = {'ADCP','RACMO','MAR','MERRA',['ICE (' upper(userdata) ')']};
end

% FIGURE 1 = RUNOFF (options: 'raw','mean','members','sensitivity','surf')
if opts.simmodel == "skinmodel"
   h1 = plotRunoff( ...
      Runoff,Discharge,Catchment, ...
      'plotsurf',true, ...
      'legendtext',ltext, ...
      'sitename',sitename, ...
      'refstart',true ...
      );
else
   if ismember(sitename, {'slv1','slv2'})

      h1 = plotRunoff( ...
         Runoff,Discharge,Catchment, ...
         'plotensemble',false, ...
         'legendtext',ltext, ...
         'sitename',sitename, ...
         'refstart',true ...
         );

      h1 = plotRunoff( ...
         Runoff,Discharge,Catchment, ...
         'plotensemble',false, ...
         'legendtext',ltext, ...
         'sitename',sitename, ...
         'refstart',false ...
         );
   else
      h1 = plotRunoff( ...
         Runoff,Discharge,Catchment, ...
         'plotensemble',false, ...
         'legendtext',ltext, ...
         'sitename',sitename ...
         );
   end
end

% PLOT ABLATION
plotPromice(AblationHourly, 'refstart', ...
   datetime(str2double(simyears), 6, 1, 0, 0, 0, 'TimeZone', 'UTC'));
plotPromice(AblationDaily, 'refstart', ...
   datetime(str2double(simyears), 6, 1, 0, 0, 0, 'TimeZone', 'UTC'));

plotPromice(AblationHourly, 'refstart', ...
   datetime(str2double(simyears), 7, 1, 0, 0, 0, 'TimeZone', 'UTC'));
plotPromice(AblationDaily, 'refstart', ...
   datetime(str2double(simyears), 7, 1, 0, 0, 0, 'TimeZone', 'UTC'));

% plot enbal - need option to plot met station data when forcing is rcm
plotEnbal(ice1, met);

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
