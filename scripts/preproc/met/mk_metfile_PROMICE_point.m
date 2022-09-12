clean

% these are finished running for KAN-M and KAN-L w/wo MODIS albedo 2009-18

saveData    = false;
% yrs         = 2009:2018;
% siteName    = 'KAN_L';
yrs         = 2016;
siteName    = 'KAN_M';
timeScale   = 'hourly';
numYears    = length(yrs);
plotMet     = false;            % plot the data?
tinterp     = true;
newdt       = '15m';

% fill_albedo = fill missing values. fill_winter means set NDJJM to 0.8. if
% fill_winter == false, then the first / last valid albedo value is used to
% backfill the missing witner values, but sometimes the first / last valid
% value is very low (summer value), so default to fill_winter = true 
fillAlbedo  = true;
fillWinter  = true;
plotFill    = false;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% set paths and load data
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p.data  = '/Users/coop558/mydata/geus/aws/v3/';
p.save  = setpath('GREENLAND/icemodel/input/met/');
p.modis = '/Users/coop558/mydata/geus/albedo/raw/'; % for swapping out albedo

% Read the PROMICE data. See notes for using pre-saved PROMICE .mat file.
[~,met] = readPromiceAWS(p.data,siteName,timeScale);

% i pull out the structure components b/c I use 'met' for the metfile
data    = met.data;
ind     = met.inds;
t       = met.t;
header  = string(met.header); clear met;
Tmet    = datetime(t,'ConvertFrom','datenum');

load('projsipsn.mat');

% these header names don't match my conventions, and it's easier to
% manually fix here than add if/else statments below
ind.swd     = ind.swd_corrected;
ind.tsfc    = ind.tsrf;

% more compact for filenames
siteName2   = strrep(siteName,'_','');

% albedo requires special pre-processing
if fillAlbedo == true
    data = fillPromiceAlbedo(data,Tmet,header,fillWinter,plotFill);
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% pull out the station data for the period of interest
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% use one of the mar metfiles as a template
load([p.save 'met_behar_MAR_2016_1hr.mat']);
vars    =   met.Properties.VariableNames;
units   =   met.Properties.VariableUnits; clear met

for n = 1:numYears

    yy      = yrs(n);
    yys     = int2str(yy);
    T       = datetime(yy,1,1,0,0,0):hours(1):datetime(yy,12,31,23,0,0);
    Tday    = datetime(yy,1,1):days(1):datetime(yy,12,31);
    timeIdx = ismember(Tmet,T);
    
    % set file names that change each year
    fmodis  = [p.modis 'Greenland_Reflectivity_' yys '_5km_C6.nc'];
    fsave   = ['met_' siteName2 '_' siteName2 '_'];
    
    % calculate positional information
    lat     = nanmean(data(timeIdx,ind.lat));
    lon     = nanmean(data(timeIdx,ind.lon)); 
    [x,y]   = projfwd(projsipsn,lat,-lon);
    elev    = nanmean(data(timeIdx,ind.elev));
    slope   = nan;

   for m = 1:numel(vars)
      
      thisVar     =  vars{m};
      
      % COPY THE DATA
      if isfield(ind,thisVar)
         varIdx         = ind.(thisVar);
         met.(thisVar)  = data(timeIdx,varIdx);
      else
         met.(thisVar)  = nan(size(met.tair));
      end
         
      % UNIT CONVERSIONS
      if thisVar == "tair" || thisVar == "tsfc"
         met.(thisVar)  = met.(thisVar) + 273.15;        % C -> K
      elseif thisVar == "psfc"
         met.(thisVar)  = met.(thisVar) .*100;           % hPa -> Pa
      end
      
      % the ones set nan are: snow, rain, melt, runoff, snowd, smb
   end
    
   
   met.date    = tocolumn(datenum(T));
    
% add modis albedo
   LON         = double(ncread(fmodis,'lon'));
   LAT         = double(ncread(fmodis,'lat'));
   [X,Y]       = projfwd(projsipsn,LAT,LON);
   [xs,xc]     = ncrowcol(X,X,Y,x,y);      % x/y start for ncread
   vs          = [xs 1];                   % var start for ncread
   vc          = [xc numel(T)/24];         % var count
   alb         = squeeze(double(ncread(fmodis,'albedo',vs,vc)));
   met.MODIS   = tocolumn(interp1(Tday,alb,T));
   
% perform some basic checks    
    met         = metchecks(met,false);  % false = don't plot

% convert to timetable    
    met     = struct2table(met);
    met     = table2timetable(met,'RowTimes',T);

% define units for the scalar quantities    
    scalarunits = ["m" "m" "oN" "oW" "m asl" "m/m"];
    
% add scalar properties
    met     =   addprop(met, {'X','Y','Lat','Lon','Elev','Slope',       ...
                              'ScalarUnits'},{'table','table',          ...
                              'table','table','table','table','table'});

% update the properties
%   met.Properties.VariableNames                = vars;
    met.Properties.VariableUnits                = units;
    met.Properties.CustomProperties.X           = x;
    met.Properties.CustomProperties.Y           = y;
    met.Properties.CustomProperties.Lat         = lat;
    met.Properties.CustomProperties.Lon         = lon;
    met.Properties.CustomProperties.Elev        = elev;
    met.Properties.CustomProperties.Slope       = slope;
    met.Properties.CustomProperties.ScalarUnits = scalarunits; 

    
% % % 
% metday = retime(met,'daily','mean');
% figure; plot(Tday,alb); hold on; plot(metday.Time,metday.albedo);
% legend('modis','kanm')
% 
% test = metday;
% test.modis = alb;
% 
% test.albedo = smoothdata(test.albedo,'movmean',30);
% test.modis = smoothdata(test.modis,'movmean',30);
% 
% figure; plot(test.Time,test.modis); hold on; plot(test.Time,test.albedo);
% legend('modis','kanm')
% 
% figure; plot(met.Time,met.albedo);
% 
% alb = round(alb,3);
% idx = find(alb==0.695);
% 
% plot(Tday(idx),alb(idx))

% % % 

    if saveData == true
        save([p.save fsave yys '_1hr.mat'],'met')

    % interpolate, re-check out-of-bound values, and save
        if tinterp == true
            met      = interpMet(met,newdt);
            met.date = datenum(met.Time);
            met      = metchecks(met,false); % re-check out-of-bound values

            save([p.save fsave yys '_' newdt '.mat'],'met');
        end

    end
    
    clear met
end

% % % % % % % % % % % % % % % % % % % % % % % 
% these header names (on ind.) are from the .nc file in the
% read_save_PROMICE script that used the .nc files just to get headers.
% with the read_promice_aws script, i use my custom header renaming

%     met.tair    = data(idx,ind.air_temperature) + 273.16;       % C -> K
%     met.swd     = data(idx,ind.shortwave_radiation_down_cor);
%     met.lwd     = data(idx,ind.longwave_radiation_down);
%     met.albedo  = data(idx,ind.albedo_theta);
%     met.snow    = nan(size(met.albedo));
%     met.rain    = nan(size(met.albedo));
%     met.melt    = nan(size(met.albedo));
%     met.runoff  = nan(size(met.albedo));
%     met.shf     = data(idx,ind.sensible_heat_flux);
%     met.lhf     = data(idx,ind.latent_heat_flux);
%     met.wspd    = data(idx,ind.wind_speed);
%     met.wdir    = data(idx,ind.wind_direction);
%     met.psfc    = data(idx,ind.air_pressure).*100;              % hPa -> Pa
%     met.rh      = data(idx,ind.relative_humidity);
%     met.snowd   = nan(size(met.albedo));
%     met.tsfc    = data(idx,ind.surface_temp)  + 273.16;         % C -> K
%     met.cfrac   = data(idx,ind.cloudcover);
%     met.date    = tocolumn(datenum(T));


% % % % % % % % % % % % % % % 



% % % % % % % % % 
% originally I read in the KAN_L and KAN_M pre-saved .mat files. Then I
% switched out so I could read in PROMICE and pick station. Below shows how
% that works, and after that shows they are identcial. But then I switched
% to the function, because it reads teh raw .txt files and is better
% practice overall.

% % this loads the pre-saved PROMICE data structure
% p.data  = setpath('GREENLAND/field/meteo/PROMICE/data/2018/');
% load([p.data 'PROMICE_2008_2018_' timescale '.mat']);
% imet    = contains(stations,site);
% data    = data{imet};
% t       = t{imet}; % clear header stations imet;

% % this shows that the processed kan-l and new extract are identical,
% where the new extract is using PROMICE_ data structure with all stations
% figure; 
% for n = 1:size(data,2)-6
% 
%     if strcmp(site,'KAN_L')
%         plot(met.data(:,n+6),data(:,n+6)); addOnetoOne(':');
%         title(met.header{n+6});
%     elseif strcmp(site,'KAN_M')
%         plot(kanm.data(:,n+6),data(:,n+6)); addOnetoOne(':');
%         title(kanm.header{n+6});
%     end
%     pause
%     clf
% end

% when using the above version, I needed these lines immediately after
% loading the data. I also had a switch statement to read kanl or kanm
% data    = met.data;
% ind     = met.inds;
% t       = met.t;    clear met
% Tmet    = datetime(t(:,1),t(:,2),t(:,3),t(:,4),t(:,5),0);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% the old way
% met = [Y,M,D,H,tair,rh,wspd,wdir,ppt,qsi,qli,albedo,press];
% ice melt input is: [yyyy mm dd hh tair rh wspd wdir pppt]
% mk_icemelt_station_file(  yyi,mmi,ddi,hhi,dt,yyf,mmf,ddf,hhf, ...
%                           tair,rh,wspd,wdir,ppt,fname);

    
%     % start time
%     yyi         =   yrs(n);
%     mmi         =   1;
%     ddi         =   1;
%     hhi         =   0;
%     mni         =   0;
%     sci         =   0;
% 
%     % end time
%     yyf         =   yrs(n);
%     mmf         =   12;
%     ddf         =   31;
%     hhf         =   23;
%     mnf         =   59;
%     scf         =   60;
% 
%     % timestep, in hours
%     dt          =   1/24;
% 
%     % build calendar
%     dni         =   datenum(yyi,mmi,ddi,hhi,mni,sci);
%     dnf         =   datenum(yyf,mmf,ddf,hhf,mnf,scf);
%     tmodel      =   dni:dt:dnf;
% 
%     % get the si, ei on the Kan-M time vector
%     [si,ei]     =   dateInds(yyi,mmi,ddi,hhi,yyf,mmf,ddf,hhf,t);
% 
%     % % below are items needed as input to icemelt.f, I will pick up on creating
%     % % this file once the promice data is obtained
%     % 
%     % elev        =   nanmean(data(idx,ind.Elev));
%     % lat         =   nanmean(data(idx,ind.Lat))/100;
%     % tair        =   nanmean(data(idx,ind.Tair));
%     % tice1       =   nanmean(data(idx,ind.ice_temp_01));
%     % tice2       =   nanmean(data(idx,ind.ice_temp_02));
%     % tice3       =   nanmean(data(idx,ind.ice_temp_03));
%     % tice4       =   nanmean(data(idx,ind.ice_temp_04));
%     % tice5       =   nanmean(data(idx,ind.ice_temp_05));
%     % tice6       =   nanmean(data(idx,ind.ice_temp_06));
%     % tice7       =   nanmean(data(idx,ind.ice_temp_07));
%     % tice8       =   nanmean(data(idx,ind.ice_temp_08));
% 
%     %% pull out the data needed for the met file
%     id          =   1;
%     lat         =   nanmean(data(idx,ind.latitude_GPS));
%     lon         =   nanmean(data(idx,ind.longitude_GPS)); 
%     lon         =   -1*lon;
%     [x,y]       =   projfwd(proj_nps,lat,lon);
%     elev        =   nanmean(data(idx,ind.elevation));
%     tair        =   data(idx,ind.air_temperature);
%     rh          =   data(idx,ind.relative_humidity);
%     wspd        =   data(idx,ind.wind_speed);
%     wdir        =   data(idx,ind.wind_direction);
%     ppt         =   zeros(size(wdir));
%     qsi         =   data(idx,ind.shortwave_radiation_down_cor);
%     qli         =   data(idx,ind.longwave_radiation_down);
%     albedo      =   data(idx,ind.albedo_theta);
%     press       =   data(idx,ind.air_pressure);
%     % tsfc        =   data(idx,ind.surface_temp);
% 
% %     fname       =   [p.save 'met.dat'];




% Tpor    = datetime(yrs(1),1,1,0,0,0):hours(1):datetime(yrs(end),12,31,23,0,0);
% idx     = ismember(Tmet,Tpor);
% alb     = data(idx,ind.albedo_theta);
% si      = 1:1000;
% ei      = length(alb)-1000:length(alb);
% alb(si) = 0.8;
% alb(ei) = 0.8;
% albq    = fillmissing(alb,'linear');
% 
% % figure; plot(Tpor,albq); hold on; plot(Tpor,alb);