clean

save_data           =   1;
fname               =   'met.dat';
%%

homepath            =   pwd;

if strcmp(homepath(2:6),'Users')
    path.data       =   ['/Users/mattcooper/Dropbox/CODE/MATLAB/GREENLAND/' ...
                            'field/meteo/PROMICE/data/2018/'];
    path.save       =   ['/Users/mattcooper/Dropbox/CODE/MATLAB/GREENLAND/' ...
                            'runoff/icemodel/data/preprocess/output/2015_2016/'];
elseif strcmp(homepath(10:16),'mcooper')    
end

%%
load([path.data 'kan_M_2008_2018_hourly.mat']);
load('proj_nps');

%% pull out the kan m data for the period of interest

% start time
yyi         =   2015;
mmi         =   1;
ddi         =   1;
hhi         =   0;
mni         =   0;
sci         =   0;

% end time
yyf         =   2016;
mmf         =   12;
ddf         =   31;
hhf         =   23;
mnf         =   59;
scf         =   60;

% timestep, in hours
dt          =   1/24;

% build calendar
dni         =   datenum(yyi,mmi,ddi,hhi,mni,sci);
dnf         =   datenum(yyf,mmf,ddf,hhf,mnf,scf);
tmodel      =   dni:dt:dnf;

% get the si, ei on the Kan-M time vector
[si,ei]     =   dateInds(yyi,mmi,ddi,hhi,yyf,mmf,ddf,hhf,t);

% % below are items needed as input to icemelt.f, I will pick up on creating
% % this file once the promice data is obtained
% 
% elev        =   nanmean(data(si:ei,ind.Elev));
% lat         =   nanmean(data(si:ei,ind.Lat))/100;
% tair        =   nanmean(data(si:ei,ind.Tair));
% tice1       =   nanmean(data(si:ei,ind.ice_temp_01));
% tice2       =   nanmean(data(si:ei,ind.ice_temp_02));
% tice3       =   nanmean(data(si:ei,ind.ice_temp_03));
% tice4       =   nanmean(data(si:ei,ind.ice_temp_04));
% tice5       =   nanmean(data(si:ei,ind.ice_temp_05));
% tice6       =   nanmean(data(si:ei,ind.ice_temp_06));
% tice7       =   nanmean(data(si:ei,ind.ice_temp_07));
% tice8       =   nanmean(data(si:ei,ind.ice_temp_08));

%% pull out the data needed for the met file
id          =   1;
lat         =   nanmean(data(si:ei,ind.latitude_GPS));
lon         =   nanmean(data(si:ei,ind.longitude_GPS)); 
lon         =   -1*lon;
[x,y]       =   projfwd(proj_nps,lat,lon);
elev        =   nanmean(data(si:ei,ind.elevation));
tair        =   data(si:ei,ind.air_temperature);
rh          =   data(si:ei,ind.relative_humidity);
wspd        =   data(si:ei,ind.wind_speed);
wdir        =   data(si:ei,ind.wind_direction);
ppt         =   zeros(size(wdir));

%% make the micromet file

% reset dt to 1 for 1 hour
dt          =   1;

% ice melt input is:
% [yyyy mm dd hh tair rh wspd wdir pppt]

% note - icemelt station file has fewer requirements than micromet
fname               =   [path.save fname];

if save_data == 1
    mk_icemelt_station_file(    yyi,mmi,ddi,hhi,dt,yyf,mmf,ddf,hhf, ...
                                tair,rh,wspd,wdir,ppt,fname);
end












