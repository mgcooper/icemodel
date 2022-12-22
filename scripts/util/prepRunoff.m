function [Runoff,Discharge,Catchment] = prepRunoff(opts,ice1)

site = opts.sitename;
yyyy = num2str(opts.simyears(1));

patheval = getenv('ICEMODELDATAPATH');
pathdata = [getenv('ICEMODELINPUTPATH') 'userdata/'];

% for the purpose of prepRunoff, it might make sense to just swap the siteName
if strcmpi(site,'hills')
   site = 'behar';
   disp('site=hills, using behar runoff')
end
if strcmpi(site,'kanm')
   site = 'behar';
   disp('site=kanm, using behar runoff')
end
if strcmpi(site,'kanl')
   site = 'ak4';
   disp('site=kanl, using ak4 runoff')
end

% read in the observed discharge and catchment area. need the if/else for behar
% until the 2015/16 data is combined into one file
if strcmpi(site,'behar')
   load([patheval '/Q_behar_' yyyy '.mat'],'Q','sf');
else
   load([patheval '/Q_' site],'Q','sf');
end

Discharge = Q;

% UpperBasin has no min/max, but the uncertainty was estimated to be +/- 5%
if ~isfield(sf,'min')
   sf.min.ease.area = 0.95*sf.med.ease.area;
end
if ~isfield(sf,'max')
   sf.max.ease.area = 1.05*sf.med.ease.area;
end
Catchment = sf;
Catchment.sitename = site;

% read in the RCM model data. for years with missing data, use MAR
MAR = load([pathdata 'mar_' site '_' yyyy '.mat'],'Data');
if strcmp(site,'ak4')
   if str2double(yyyy)>=2012 && str2double(yyyy)<=2016
      RAC = load([pathdata 'racmo_' site '_' yyyy '.mat'],'Data');
      MER = load([pathdata 'merra_' site '_' yyyy '.mat'],'Data');
   else
      RAC = MAR;
      RAC.Data.runoff   = nan.*RAC.Data.runoff;
      RAC.Data.melt     = nan.*RAC.Data.melt;

      MER = MAR;
      MER.Data.runoff   = nan.*MER.Data.runoff;
      MER.Data.melt     = nan.*MER.Data.runoff;
      MER.Data.rain     = nan.*MER.Data.runoff;
   end
else
   RAC = load([pathdata 'racmo_' site '_' yyyy '.mat'],'Data');
   MER = load([pathdata 'merra_' site '_' yyyy '.mat'],'Data');
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

% check if there is leap data
if isleap(str2double(yyyy))
   T = ice1.Time;
   Tleap = datetime(year(T(1)),1,1,0,0,0):hours(1):datetime(year(T(1)),12,31,23,0,0);
   feb29 = month(Tleap)==2 & day(Tleap)==29;
   fields = fieldnames(Runoff);
   for n = 1:numel(fields)
      dat = Runoff.(fields{n});
      if size(dat,1) == numel(Tleap)
         dat = rmleapinds(dat,Tleap);
         Runoff.(fields{n}) = dat;
      end
   end
end

% convert Runoff to a timetable
Runoff.Time     = ice1.Time;
Runoff          = table2timetable(struct2table(Runoff));

% % shouldn't be needed anymore
% % if opts.dt > 3600
% %     ice1    = retime(ice1,datetime(mar.dates,'ConvertFrom','datenum'),  ...
% %                     'linear');
% % end

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