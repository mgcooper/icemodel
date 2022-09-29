function Ablation = prepAblation(opts,ice1,timescale)
   
%    addParameter(  p, 'timescale',      '',         @(x)ischar(x)        );

   site    = opts.sitename;
   yyyy    = num2str(opts.simyears(1));
   
   patheval = getenv('ICEMODELDATAPATH');
   pathdata = [getenv('ICEMODELINPUTPATH') 'userdata/'];
   
   % load the promice ablation data
   if contains(site,{'behar','slv1','slv2','kanm','KANM'})
      aws = 'kanm';
   elseif contains(site,{'upperbasin','ak4','kanl','KANL'})
      aws = 'kanl';
   end
   load([patheval aws '_ablation_' timescale '.mat'],'ablation');
   
   % pull out the ablation for the simulation time
   iplot       = isbetween(ablation.Time,ice1.Time(1),ice1.Time(end));
   Ablation    = ablation(iplot,:);
   
   % read in the RCM model data
   if exist([pathdata 'mar_' site '_' yyyy '.mat'],'file') == 2
      MAR               = load([pathdata 'mar_' site '_' yyyy '.mat'],'Data');
   else
      MAR.Data.runoff   = nan(size(Ablation.ablation));
      MAR.Data.melt     = nan(size(Ablation.ablation));
   end
   if exist([pathdata 'racmo_' site '_' yyyy '.mat'],'file') == 2
      RAC               = load([pathdata 'racmo_' site '_' yyyy '.mat'],'Data');
   else % create a dummy timetable for compatibility but set nan 
      RAC.Data          = MAR.Data;
      RAC.Data.runoff   = nan.*RAC.Data.runoff;
      RAC.Data.melt     = nan.*RAC.Data.melt;
   end
   if exist([pathdata 'merra_' site '_' yyyy '.mat'],'file') == 2
      MER               = load([pathdata 'merra_' site '_' yyyy '.mat'],'Data');
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

   
      
   % this data is no longer used
   % pathaws = setpath('runoff/data/ablation/mat/','project');
   % load([pathaws 'promice_ablation_annual.mat'],'data');
   % iyear  = find(data.Year == str2double(yyyy));
   % AnnualAblation  = data.KAN_M(iyear);
   % AnnualAblation  = data.KAN_L(iYear);
   
end
