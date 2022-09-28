function [Runoff,Discharge,Catchment] = prepRunoff(opts,ice1)
    
   site    = opts.sitename;
   yyyy    = opts.yyyy;
   
   % for the purpose of prepRunoff, it might make sense to just swap the
   % siteName
   if strcmpi(site,'hills'); 
      site = 'behar'; 
      disp('site=hills, using behar runoff')
   end
   if strcmpi(site,'kanm'); 
      site = 'behar'; 
      disp('site=kanm, using behar runoff')
   end
   if strcmpi(site,'kanl'); 
      site = 'ak4'; 
      disp('site=kanl, using ak4 runoff')
   end
%    if strcmpi(site,'kanl'); site = 'upperBasin'; end

   p.data   = '/Users/coop558/mydata/';
   p.mar    = [p.data 'mar3.11/matfiles/' site '/data/'];
   p.rac    = [p.data 'racmo2.3/matfiles/' site '/data/'];
   p.mer    = [p.data 'merra2/1hrly/matfiles/' site '/data/'];
   p.obs    = setpath(['GREENLAND/runoff/data/icemodel/eval/' site '/']);
% % this was an attemp to get rid of the siteName swapping   
%    % this finds the nearest site with runoff data
%    if strcmp(site,'behar') || strcmpi(site,'kanm')
%       p.obs = setpath('GREENLAND/runoff/data/icemodel/eval/behar/');
%       load([p.obs opts.yyyy '/Q_behar.mat'],'Q','sf');
%    elseif strcmpi(site,'kanl')
%       p.obs = setpath('GREENLAND/runoff/data/icemodel/eval/ak4/');
%       load([p.obs 'Q_ak4.mat'],'Q','sf');
%    end
   
% % if not using the commented out, need this for behar
   if strcmpi(site,'behar')
      p.obs = setpath('runoff/data/icemodel/eval/behar/','project');
      p.obs = [p.obs opts.yyyy '/'];
   end
   
% read in the observed discharge and catchment area
   load([p.obs 'Q_' site],'Q','sf');
   Discharge   = Q;

% UpperBasin has no min/max, but the uncertainty was estimated to be +/- 5% 
    if ~isfield(sf,'min')
        sf.min.ease.area = 0.95*sf.med.ease.area;
    end
    if ~isfield(sf,'max')
        sf.max.ease.area = 1.05*sf.med.ease.area;
    end
    Catchment   = sf;
    Catchment.sitename = site;

% read in the RCM model data
   MAR     = load([p.mar 'mar_' site '_' yyyy '.mat'],'Data');
   if strcmp(site,'ak4')
      if str2double(yyyy)>=2012 && str2double(yyyy)<=2016
         RAC   = load([p.rac 'racmo_' site '_' yyyy '.mat'],'Data');
         MER   = load([p.mer 'merra_' site '_' yyyy '.mat'],'Data');
      else
         RAC               = MAR;
         RAC.Data.runoff   = nan.*RAC.Data.runoff;
         RAC.Data.melt     = nan.*RAC.Data.melt;
      
         MER               = MAR;
         MER.Data.runoff   = nan.*MER.Data.runoff;
         MER.Data.melt     = nan.*MER.Data.runoff;
         MER.Data.rain     = nan.*MER.Data.runoff;
      end
   else
      RAC   = load([p.rac 'racmo_' site '_' yyyy '.mat'],'Data');
      MER   = load([p.mer 'merra_' site '_' yyyy '.mat'],'Data');
   end
    
% for the case where an ensemble of icemodel/skinmodel runs are passed in:
    numRuns = width(ice1.runoff);
    
% combine the runoff data, convert icemodel from cumulative to incremental
    Runoff.marRunoff        = MAR.Data.runoff;
    Runoff.racmoRunoff      = RAC.Data.runoff;
    Runoff.merraRunoff      = MER.Data.runoff; 
    Runoff.icemodelRunoff   = [zeros(1,numRuns);diff(ice1.runoff)];
    Runoff.skinmodelRunoff  = [zeros(1,numRuns);diff(ice1.runoff)];
    Runoff.marMelt          = MAR.Data.melt;
    Runoff.racmoMelt        = RAC.Data.melt;
    Runoff.merraMelt        = MER.Data.runoff - MER.Data.rain;
    Runoff.icemodelMelt     = [zeros(1,numRuns);diff(ice1.melt)];
    Runoff.skinmodelMelt    = [zeros(1,numRuns);diff(ice1.melt)];
    
    % convert Runoff to a timetable
    Runoff.Time     = ice1.Time;
    Runoff          = table2timetable(struct2table(Runoff));
    
% %     % make a similar table for ablation
% %     Ablation.mar        = MAR.Data.melt;
% %     Ablation.racmo      = RAC.Data.melt;
% %     Ablation.merra      = MER.Data.runoff - MER.Data.rain;
% %     Ablation.icemodel   = [0;diff(ice1.depth_melt)];
% %     
% %     % rename 'Ablation.ablation' to Ablation.promice
% %     Ablation    = renamevars(Ablation,'ablation','promice');
    
  % merra.melt      = MER.ppt-MER.smb-MER.evap-MER.rain;
    
    


% % shouldn't be needed anymore
% % if opts.dt > 3600
% %     ice1    = retime(ice1,datetime(mar.dates,'ConvertFrom','datenum'),  ...
% %                     'linear');
% % end
% 
% % plot(T,Drn)


% 
% 
% Amin    =   sf.min.ease.area;
% Amed    =   sf.med.ease.area;
% Amax    =   sf.max.ease.area;
% 
% RmarL   =   Amin.*(cumsum(MAR.RUH)-MAR.RUH(1))./1e9;
% RmarM   =   Amed.*(cumsum(MAR.RUH)-MAR.RUH(1))./1e9;
% RmarH   =   Amax.*(cumsum(MAR.RUH)-MAR.RUH(1))./1e9;
% 
% figure;
% plot(R.Time,RmarL,':'); hold on;
% plot(R.Time,RmarH,':');
% 
% data.runoff         = data.RUH;
% data.column_runoff  = data.RUH;
% 
% mar         = data;
% racmo       = data;
% merra       = data;
% 
% h1 = plotrunoff(mar,merra,racmo,ice1,Q,sf,'raw','column');


end
    