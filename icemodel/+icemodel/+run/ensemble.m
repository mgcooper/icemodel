function ensemble(kwargs)

   arguments
      kwargs.savedata = true
      kwargs.casename = [] % use 0 for the "otherwise" case
   end
   [savedata, casename] = dealout(kwargs);

   % Each icemodel/output/<version>/<sitename>/<yearname> folder should have:
   % AWS + MAR, MERRA, MODIS, NONE, RACMO
   % MAR + AWS, MERRA, MODIS, NONE, RACMO
   % Except pre-2012 (Leverett and AK4), which won't have RACMO
   % Thus to quickly confirm all sim's are done, just count as above.
   % As of Apr 2024, all sim's are done.

   % valid userdata:
   % {'racmo','mar','merra','kanm','kanl','modis','none'}

   % valid sitenames:
   % {'upperbasin','behar','slv1','slv2','ak4'}

   % valid forcings:
   % {'MAR','KANM','KANL'}

   %% Set which predefined
   switch casename

      case 0
         % Custom configuration

         % Apr 2024 - need to use this to run any merra ones that weren't possible
         % until I made the 2009-2011 files:
         % DONE skinmodel, ak4
         % icemodel, ak4
         sitename = {'ak4'};
         forcings = {'mar', 'kanl'};
         userdata = {'merra'};
         uservars = {'albedo'};
         smbmodel = {'icemodel'};
         simyears = 2009:2011;

      case 1
         % DONE upperbasin, 2016, mar, kanl, icemodel
         sitename = {'upperbasin'};
         forcings = {'mar', 'kanl'};
         userdata = {'none', 'mar', 'merra', 'modis', 'kanl', 'racmo'};
         uservars = {'albedo'};
         smbmodel = {'icemodel', 'skinmodel'};
         simyears = 2016:2016;

      case 2
         % DONE behar, 2016, mar, kanm
         sitename = {'behar'};
         forcings = {'mar', 'kanm'};
         userdata = {'none', 'mar', 'merra', 'modis', 'kanm', 'racmo'};
         uservars = {'albedo'};
         smbmodel = {'icemodel', 'skinmodel'};
         simyears = 2016:2016;

      case 3
         % DONE behar, slv1, slv2, 2015, mar, kanm, icemodel
         sitename = {'behar', 'slv1', 'slv2'};
         forcings = {'mar', 'kanm'};
         userdata = {'none', 'mar', 'merra', 'modis', 'kanm', 'racmo'};
         uservars = {'albedo'};
         smbmodel = {'icemodel', 'skinmodel'};
         simyears = 2015:2015;

      case 4
         % DONE ak4, 2009-2011, mar, kanl, icemodel
         sitename = {'ak4'};
         forcings = {'mar', 'kanl'};
         userdata = {'none', 'mar', 'merra', 'modis', 'kanl'};
         uservars = {'albedo'};
         smbmodel = {'skinmodel'}; % 'icemodel',
         simyears = 2009:2011;

      case 5
         % DONE ak4, 2012-2016, mar, kanl, icemodel, merra, racmo
         sitename = {'ak4'};
         forcings = {'mar', 'kanl'};
         userdata = {'none', 'mar', 'merra', 'modis', 'kanl', 'racmo'};
         uservars = {'albedo'};
         smbmodel = {'skinmodel'}; % 'icemodel',
         simyears = 2012:2016;

      otherwise

         error('unrecognized run combo')

   end

   %% Build the ensemble list
   ensemble = ensembleList( ...
      forcings, userdata, uservars, smbmodel, simyears, sitename);

   %% Run all combos
   for n = 1:ensemble.numcombos

      sitename = char(ensemble.allcombos.sitename(n));
      simyears = char(ensemble.allcombos.simyears(n));
      smbmodel = char(ensemble.allcombos.smbmodel(n));
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
      opts = icemodel.setopts(smbmodel, sitename, str2double(simyears), ...
         forcings, userdata, uservars, savedata);

      % display the run info
      disp([smbmodel ', ' sitename ', '  forcings ', ' userdata ', ' simyears])

      % RUN THE MODEL
      switch smbmodel
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
   if opts.smbmodel == "skinmodel"
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
end
