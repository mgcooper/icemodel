function Ablation = prepAblation(opts,ice1,timescale)
   
%    addParameter(  p, 'timescale',      '',         @(x)ischar(x)        );

   site    = opts.sitename;
   yyyy    = opts.yyyy;
   
%    % this will work if we have Data files for the site/year, but those are
%    % only for years with runoff. This means I need to make Data files for
%    % the met stations, which is needed anyway for a comprehensive ablation
%    % analysis
%    if strcmpi(site,'kanm'); site = 'behar'; end
%    if strcmpi(site,'kanl'); site = 'ak4'; end
%    %    if strcmpi(site,'kanl'); site = 'upperBasin'; end
   
   if strcmpi(site,'hills'); 
      site = 'behar'; 
      disp('site=hills, using behar runoff')
   end

   p.data  = '/Users/coop558/mydata/';
   p.mar   = [p.data 'mar3.11/matfiles/' site '/data/'];
   p.rac   = [p.data 'racmo2.3/matfiles/' site '/data/'];
   p.mer   = [p.data 'merra2/1hrly/matfiles/' site '/data/'];
   p.obs   = setpath(['GREENLAND/runoff/data/icemodel/eval/' site '/']);
   p.aws1  = setpath('GREENLAND/runoff/data/ablation/mat/');
   p.aws2  = setpath('GREENLAND/runoff/data/icemodel/eval/');
   

   if strcmp(site,'behar')
      p.obs = [p.obs opts.yyyy '/'];
   end
   
   % load the promice ablation data
   load([p.aws1 'promice_ablation_annual.mat'],'data');
   iyear  = find(data.Year == str2double(yyyy));
   switch site
      case {'behar','slv1','slv2','kanm','KANM'}
         if strcmp(timescale,'daily')
            load([p.aws2 'kanm/kanm_ablation_daily'],'ablation');
         elseif strcmp(timescale,'hourly')
            load([p.aws2 'kanm/kanm_ablation_hourly'],'ablation');
         end
         
         %AnnualAblation  = data.KAN_M(iYear);
         
      case {'upperBasin','ak4','kanl','KANL'}
         if strcmp(timescale,'daily')
            load([p.aws2 'kanl/kanl_ablation_daily'],'ablation');
         elseif strcmp(timescale,'hourly')
            load([p.aws2 'kanl/kanl_ablation_hourly'],'ablation');
         end
         %AnnualAblation  = data.KAN_L(iYear);
   end
   
   % pull out the ablation for the simulation time
   iplot       = isbetween(ablation.Time,ice1.Time(1),ice1.Time(end));
   Ablation    = ablation(iplot,:);
   
   % read in the RCM model data
   if exist([p.mar 'mar_' site '_' yyyy '.mat'],'file')
      MAR               = load([p.mar 'mar_' site '_' yyyy '.mat'],'Data');
   else
      MAR.Data.runoff   = nan(size(Ablation.ablation));
      MAR.Data.melt     = nan(size(Ablation.ablation));
   end
   if exist([p.rac 'racmo_' site '_' yyyy '.mat'],'file')
      RAC               = load([p.rac 'racmo_' site '_' yyyy '.mat'],'Data');
   else % create a dummy timetable for compatibility but set nan 
      RAC.Data          = MAR.Data;
      RAC.Data.runoff   = nan.*RAC.Data.runoff;
      RAC.Data.melt     = nan.*RAC.Data.melt;
   end
   if exist([p.mer 'merra_' site '_' yyyy '.mat'],'file')
      MER   = load([p.mer 'merra_' site '_' yyyy '.mat'],'Data');
   else
      MER.Data          = MAR.Data;
      MER.Data.runoff   = nan.*MER.Data.runoff;
      MER.Data.melt     = nan.*MER.Data.runoff;
      MER.Data.rain     = nan.*MER.Data.runoff;
   end
      
   % to use daily promice ablation:
   if strcmp(timescale,'daily')
      mar   = retime(MAR.Data,'daily','sum');
      rac   = retime(RAC.Data,'daily','sum');
      mer   = retime(MER.Data,'daily','sum');

      melt     = [0;diff(ice1.melt)];
      runoff   = [0;diff(ice1.runoff)];
      Time     = ice1.Time;
      ice1     = timetable(Time,melt,runoff);
      ice1     = retime(ice1,'daily','sum');

      Ablation.marRunoff        = mar.runoff;
      Ablation.racmoRunoff      = rac.runoff;
      Ablation.merraRunoff      = mer.runoff;
      Ablation.icemodelRunoff   = ice1.runoff;
      Ablation.skinmodelRunoff  = ice1.runoff;

      Ablation.marMelt          = mar.melt;
      Ablation.racmoMelt        = rac.melt;
      Ablation.merraMelt        = mer.runoff - mer.rain;
      Ablation.icemodelMelt     = ice1.melt;
      Ablation.skinmodelMelt    = ice1.melt;
      
   elseif strcmp(timescale,'hourly')
   
   % to use hourly promice ablation:   
%    % combine the runoff data, convert icemodel from cumulative to incremental
      Ablation.marRunoff        = MAR.Data.runoff;
      Ablation.racmoRunoff      = RAC.Data.runoff;
      Ablation.merraRunoff      = MER.Data.runoff;
      Ablation.icemodelRunoff   = [0;diff(ice1.runoff)];
      Ablation.skinmodelRunoff  = [0;diff(ice1.runoff)];

      Ablation.marMelt          = MAR.Data.melt;
      Ablation.racmoMelt        = RAC.Data.melt;
      Ablation.merraMelt        = MER.Data.runoff - MER.Data.rain;
      Ablation.icemodelMelt     = [0;diff(ice1.melt)];
      Ablation.skinmodelMelt    = [0;diff(ice1.melt)];
   end
   
   % rename 'Ablation.ablation' to Ablation.promice
   Ablation    = renamevars(Ablation,'ablation','promice');
   
   % merra.melt      = MER.ppt-MER.smb-MER.evap-MER.rain;
   
end
