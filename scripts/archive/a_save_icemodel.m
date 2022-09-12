clean

% this script reads in the enbal, ice1, and ice2 icemodel output

%%
save_data           =   1;

%%
homepath            =   pwd;

if strcmp(homepath(2:6),'Users')
    path.data       =   ['/Users/mattcooper/Dropbox/CODE/MATLAB/' ...
                            'GREENLAND/runoff/icemodel/model/v9/output/'];
    path.save       =   ['/Users/mattcooper/Dropbox/CODE/MATLAB/' ...
                            'GREENLAND/runoff/icemodel/model/v9/output/'];
elseif strcmp(homepath(10:16),'mcooper')    
end

%%
load([path.data 'enbal.mat']);
load([path.data 'ice1.mat']);
load([path.data 'ice2.mat']);
% load([path.data 'ice3.mat']);

%% build a calendar
dates               =   enbal(:,2:5);
yyyy                =   dates(:,1);
mm                  =   dates(:,2);
dd                  =   dates(:,3);
hh                  =   dates(:,4);
mn                  =   zeros(size(hh));
ss                  =   zeros(size(hh));
dates(:,end)        =   datenum(yyyy,mm,dd,hh,mn,ss);

% get number of days
if isleap(yyyy(1))
    ndays           =   366;
else
    ndays           =   365;
end

%% ice2

% ice 2 is nlayers(500) x nhours (24*ndays) x nvars
% the z-grid is 15 m deep with dz = 0.03 m

% reshape ice2
nlayers             =   size(ice2,1);
nhours              =   size(ice2,2);
nvars               =   size(ice2,3);
                
% pull out the vars
vars                =   {'tice','gamma','Cp_snow','h_melt','h_freeze', ...
                            'h_drain','ro_snow'};
varnames            =   {'ice temperature','thermal conductivity',      ...
                        'Specific Heat Capacity','Cumulative Melt Production', ...
                        'Cumulative Refreezing','Cumulative Drainage',  ...
                        'Ice (Snow) Density'};
units               =   {'[^oC]','[W m^{-1} K]','[J kg^{-1} K^{-1}',    ...
                        '[m]','[m]','[m]','[kg m^{-3}]'};
agg_method          =   [1;1;1;2;2;2;1]; 

for n = 1:length(vars)
    var_n           =   vars{n};
    data_n          =   ice2(:,:,n);
    data_n          =   reshape(data_n,nlayers,24,ndays);
    switch agg_method(n) % average (1) or sum (2) the hourly data?
        case 1
            dat.(var_n)     =   squeeze(mean(data_n,2));
        case 2
            dat.(var_n)     =   squeeze(sum(data_n,2));
    end
end

% put back into ice2
ice2                =   dat;
ice2.vars           =   vars;
ice2.varnames       =   varnames;
ice2.units          =   units;

%% combine the data I want to keep
            
% Save the ice data.                                
data                =   [enbal(:,6:end) ice1(:,3:7)];

                        
vars                =   {'tair','tsfc','rh','wspd','qsi','qli','qle','qh', ...
                            'qe','qc','qm','balance','albedo','qsip', ...
                            'meltdepth','meltflux','freezedepth','freezeflux'};
                        
varnames            =   {'air temp','sfc temp','rel hum','wind speed',  ...
                            'solar down','longwave down','longwave up', ...
                            'sensible flux','latent flux','conductive flux', ...
                            'melt energy','balance','albedo','solar penetration', ...
                            'melt depth','melt flux','freeze depth','freeze flux'};

units               =   {'^oC','^oC','%','m/s','W m^{-2}','W m^{-2}', ...
                            'W m^{-2}','W m^{-2}','W m^{-2}','W m^{-2}', ...
                            'W m^{-2}','W m^{-2}','-','W m^{-2}','m','m', ...
                            'm','m'};

agg_method          =   [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]; %1=avg,2=sum  

%% pull out the hourly data
for n = 1:length(vars)
    var_n           =   vars{n};
    hourly.(var_n)  =   data(:,n);
end

hourly.dates        =   dates;

%% compute daily averages/sums
                    
for n = 1:length(vars)
    var_n           =   vars{n};
    data_n          =   hourly.(var_n);
    
    % average or sum the hourly data?
    switch agg_method(n)
        case 1
            daily.(var_n)   =   nanmean(reshape(data_n,24,ndays),1);
        case 2
            daily.(var_n)   =   nansum(reshape(data_n,24,ndays),1);
    end
end

% build a calendar 
daily.dates         =   hourly.dates(1:24:end,:);

%% save the data

hourly.vars         =   vars;
hourly.varnames     =   varnames;
hourly.units        =   units;

daily.vars          =   vars;
daily.varnames      =   varnames;
daily.units         =   units;
daily.aggmethod     =   agg_method; 
daily.readme        =   'aggmethod 1 = average, 2 = sum';

icemodel.daily      =   daily;
icemodel.hourly     =   hourly;

%%
if save_data == 1
    save([path.save 'icemodel.mat'],'icemodel','ice2');
end

%%

% Enbal columns are:
% iter, Tair, Tsfc, rh, wspd, qsi, qli, qle, qh, qe, qc, qm, balance, albedo
% ice1 columns are:
% tsfc, qsi, qsip, water_depth, water_flux
% combine the variables I want
% ice2 is depth-distribution with 500 layers * 0.03 m/layer = 15 m depth

%   88  format (i6,f8.1,f10.4,f8.1,f10.4,7f9.3,f15.9,f8.2)
  
% p.iter='%i';p.tair='%9.1f';p.tsfc='%10.4f';p.rh='%9.1f';p.wspd='%10.4f';
% p.qsi='%9.3f';p.qli='%9.3f';p.qle='%9.3f';p.qh='%9.3f';p.qe='%9.3f';
% p.qc='%9.3f';p.qm='%9.3f';p.bal='%15.9f';p.alb='%8.2f\n';
% 
% fspec               =   [p.iter p.tair p.tsfc p.rh p.wspd p.qsi p.qli ...
%                             p.qle p.qh p.qe p.qc p.qm p.bal p.alb];
%                         
% fid                 =   fopen([path.data 'enbal.dat'],'r');
% test                =   fscanf(fid,fspec);

% %% if you're saving the subsnowmodel full simulation grad files
% 
% save_enbal      =       0;
% save_micromet   =       0;
% save_snowpack   =       0;
% save_snowtran   =       0;
% 
% % if your using individual variable grad files
% 
% save_albedo     =       1;      % albedo
% save_ebal       =       1;      % energy balance (should equal qm or qf?)
% save_glacmelt   =       1;      % incremental glacier melt water depth
% save_prec       =       1;      % total precipitation (liquid + solid)
% save_qc         =       1;      % conductive heat flux
% save_qe         =       1;      % latent heat flux
% save_qf         =       1;      % refreeze energy flux
% save_qh         =       1;      % sensible heat flux
% save_qle        =       1;      % outgoing longwave
% save_qli        =       1;      % incoming longwave      
% save_qm         =       1;      % melt energy flux
% save_qsi        =       1;      % incoming shortwave
% save_rain       =       1;      % rain
% save_rh         =       1;      % relative humidity
% save_ro_nsnow   =       1;      % new snow density
% save_rosnow     =       1;      % snow density
% save_runoff     =       1;      % runoff
% save_snowd      =       1;      % snow depth
% save_sprec      =       1;      % snowfall depth
% save_swe        =       1;      % swe
% save_swemelt    =       1;      % swe melt
% save_tair       =       1;      % average air temp
% save_tmin       =       0;      % min air temp
% save_tmax       =       0;      % max air temp
% save_tsfc       =       1;      % surface skin temperature
% save_wbal       =       1;      % summed water balance
% save_wdir       =       1;      % wind direction
% save_wspd       =       1;      % wind speed
% 
% 
% %% define the variables you want to delete. this means the .gdat file gets
% % deleted after read_grads reads in the data and the data is saved as .mat
% 
% delete_enbal      =       0;
% delete_micromet   =       0;
% delete_snowpack   =       0;
% delete_snowtran   =       0;
% 
% delete_albedo     =       0;      % albedo
% delete_ebal       =       0;      % energy balance (should equal qm or qf?)
% delete_glacmelt   =       0;      % incremental glacier melt water depth
% delete_prec       =       0;      % total precipitation (liquid + solid)
% delete_qc         =       0;      % conductive heat flux
% delete_qe         =       0;      % latent heat flux
% delete_qf         =       0;      % refreeze energy flux
% delete_qh         =       0;      % sensible heat flux
% delete_qle        =       0;      % outgoing longwave
% delete_qli        =       0;      % incoming longwave      
% delete_qm         =       0;      % melt energy flux
% delete_qsi        =       0;      % incoming shortwave
% delete_rain       =       0;      % rain
% delete_rh         =       0;      % relative humidity
% delete_ro_nsnow   =       0;      % new snow density
% delete_rosnow     =       0;      % snow density
% delete_runoff     =       0;      % runoff
% delete_snowd      =       0;      % snow depth
% delete_sprec      =       0;      % snowfall depth
% delete_swe        =       0;      % swe
% delete_swemelt    =       0;      % swe melt
% delete_tair       =       0;      % average air temp
% delete_tmin       =       0;      % min air temp
% delete_tmax       =       0;      % max air temp
% delete_tsfc       =       0;      % surface skin temperature
% delete_wbal       =       0;      % summed water balance
% delete_wdir       =       0;      % wind direction
% delete_wspd       =       0;      % wind speed
% %% Round the data
% 
% % It's important to round the data to the least pptision desired to reduce
% % file size. Here you set the pptision in units of 10^n i.e. 10^-1 is 0.1
% % and 10^3 is 1000, so it gets rounded to that placeholder
% 
% round_albedo     =       -2;      % albedo
% round_ebal       =       -4;      % energy balance (should equal qm or qf?)
% round_glacmelt   =       -3;      % incremental glacier melt water depth
% round_prec       =       -3;      % total precipitation (liquid + solid)
% round_qc         =       -1;      % conductive heat flux
% round_qe         =       -1;      % latent heat flux
% round_qf         =       -1;      % refreeze energy flux
% round_qh         =       -1;      % sensible heat flux
% round_qle        =       -1;      % outgoing longwave
% round_qli        =       -1;      % incoming longwave      
% round_qm         =       -1;      % melt energy flux
% round_qsi        =       -1;      % incoming shortwave
% round_rain       =       -3;      % rain
% round_rh         =       0;      % relative humidity
% round_ro_nsnow   =       0;      % new snow density
% round_rosnow     =       0;      % snow density
% round_runoff     =       -3;      % runoff
% round_snowd      =       -3;      % snow depth
% round_sprec      =       -3;      % snowfall depth
% round_swe        =       -3;      % swe
% round_swemelt    =       -3;      % swe melt
% round_tair       =       -1;      % average air temp
% round_tmin       =       -1;      % min air temp
% round_tmax       =       -1;      % max air temp
% round_tsfc       =       -1;      % surface skin temperature
% round_wbal       =       -4;      % summed water balance
% round_wdir       =       0;      % wind direction
% round_wspd       =       -1;      % wind speed
% 
% %%
% ndatasets       =   length(datafolds);
% for n = 1:ndatasets
%     
%     datafold    =   datafolds{n};
%     savefold    =   savefolds{n};
%     cd(datafold)
%     
%     if ~exist([datafold 'enbal.ctl'],'file')
%         copyfile([ctlfold 'enbal.ctl']);
%     end
%     if ~exist([datafold 'micromet.ctl'],'file')
%         copyfile([ctlfold 'micromet.ctl']);
%     end
%     if ~exist([datafold 'snowpack.ctl'],'file')
%         copyfile([ctlfold 'snowpack.ctl']);
%     end
%     if ~exist([datafold 'snowtran.ctl'],'file')
%         copyfile([ctlfold 'snowtran.ctl']);
%     end
%     
%     % MICROMET VARIABLES
%     %     ta          1  0 air temperature (deg C)
%     %     rh          1  0 relative humidity (%)
%     %     u           1  0 meridional wind component (m/s)
%     %     v           1  0 zonal wind component (m/s)
%     %     wspd        1  0 wind speed (m/s)
%     %     wdir        1  0 wind direction (0-360, true N)
%     %     qsi         1  0 incoming solar radiation reaching the surface (W/m2)
%     %     qli         1  0 incoming longwave radiation reaching the surface (W/m2)
%     %     prec        1  0 precipitation (m/time_step)
%     
%     if save_micromet == 1
%         
%         temp                =   squeeze(read_grads('micromet.ctl','ta'));
%         micromet.ta         = 	squeeze(temp(row,col,timestep)); clear temp
%         temp                =   squeeze(read_grads('micromet.ctl','rh'));
%         micromet.rh         = 	squeeze(temp(row,col,timestep)); clear temp
%         temp                = 	squeeze(read_grads('micromet.ctl','u'));
%         micromet.uwind  	=   squeeze(temp(row,col,timestep)); clear temp
%         temp                = 	squeeze(read_grads('micromet.ctl','v'));
%         micromet.vwind      =   squeeze(temp(row,col,timestep)); clear temp
%         temp                = 	squeeze(read_grads('micromet.ctl','wspd'));
%         micromet.wspd       =   squeeze(temp(row,col,timestep)); clear temp
%         temp                = 	squeeze(read_grads('micromet.ctl','wdir'));
%         micromet.wdir 		=   squeeze(temp(row,col,timestep)); clear temp
%         temp                = 	squeeze(read_grads('micromet.ctl','qsi'));
%         micromet.qsi        =   squeeze(temp(row,col,timestep)); clear temp
%         temp                = 	squeeze(read_grads('micromet.ctl','qli'));
%         micromet.qli        =   squeeze(temp(row,col,timestep)); clear temp
%         temp                = 	squeeze(read_grads('micromet.ctl','prec'));
%         micromet.prec       =   squeeze(temp(row,col,timestep)); clear temp
%     end
%     
%     if delete_micromet == 1
%         delete('micromet.gdat');
%     end
%         
%     % ENBAL VARIABLES
%     %     ta          1  0 air temperature (deg C)
%     %     tsfc        1  0 surface (skin) temperature (deg C)
%     %     qsi         1  0 incoming solar radiation reaching the surface (W/m2)
%     %     qli         1  0 incoming longwave radiation reaching the surface (W/m2)
%     %     qle         1  0 emitted longwave radiation (W/m2)
%     %     qh          1  0 sensible heat flux (W/m2)
%     %     qe          1  0 latent heat flux (W/m2)
%     %     qc          1  0 conductive heat flux (W/m2)
%     %     qm          1  0 melt energy flux (W/m2)
%     %     albedo      1  0 albedo (0-1)
%     %     ebal        1  0 energy balance error (W/m2)
%     
%     if save_enbal == 1
%         temp 				=   squeeze(read_grads('enbal.ctl','ta'));
%         enbal.ta            = 	squeeze(temp(row,col,timestep)); clear temp
%         temp            	=   squeeze(read_grads('enbal.ctl','tsfc'));
%         enbal.tsfc          = 	squeeze(temp(row,col,timestep)); clear temp
%         temp                = 	squeeze(read_grads('enbal.ctl','qsi'));
%         enbal.qsi           =   squeeze(temp(row,col,timestep)); clear temp
%         temp                = 	squeeze(read_grads('enbal.ctl','qli'));
%         enbal.qli           =   squeeze(temp(row,col,timestep)); clear temp
%         temp                = 	squeeze(read_grads('enbal.ctl','qle'));
%         enbal.qle           =   squeeze(temp(row,col,timestep)); clear temp
%         temp            	= 	squeeze(read_grads('enbal.ctl','qh'));
%         enbal.qh            =   squeeze(temp(row,col,timestep)); clear temp
%         temp             	= 	squeeze(read_grads('enbal.ctl','qe'));
%         enbal.qe            =   squeeze(temp(row,col,timestep)); clear temp
%         temp                = 	squeeze(read_grads('enbal.ctl','qc'));
%         enbal.qc            =   squeeze(temp(row,col,timestep)); clear temp
%         temp                = 	squeeze(read_grads('enbal.ctl','qm'));
%         enbal.qm            =   squeeze(temp(row,col,timestep)); clear temp
%         temp                = 	squeeze(read_grads('enbal.ctl','albedo'));
%         enbal.albedo        =   squeeze(temp(row,col,timestep)); clear temp
%         temp                = 	squeeze(read_grads('enbal.ctl','ebal'));
%         enbal.error         =   squeeze(temp(row,col,timestep)); clear temp
%         % derived variables
%         enbal.qsr           =   -1.*(1-enbal.albedo) .* enbal.qsi;
%         enbal.qsnet         =   enbal.qsi + enbal.qsr;
%         enbal.qlnet         =   enbal.qli + enbal.qle;
%         enbal.qnet          =   enbal.qsnet + enbal.qlnet;
%         % compute the energy balance
%         enbal.ebal          =   enbal.qnet + enbal.qe + enbal.qh + enbal.qc;
%     end
%     
%     if delete_enbal == 1
%         delete('enbal.gdat');
%     end
%         
%     % SNOWPACK VARIABLES
%     %     snowd       1  0 snow depth (m)
%     %     rosnow      1  0 snow density (kg/m3)
%     %     swed        1  0 snow-water-equivalent depth (m)
%     %     runoff      1  0 runoff from base of snowpack (m/time_step)
%     %     rain        1  0 liquid precipitation (m/time_step)
%     %     sprec       1  0 solid precipitation (m/time_step)
%     %     qcs         1  0 canopy sublimation (m/time_step)
%     %     canopy      1  0 canopy interception store (m)
%     %     sumqcs      1  0 summed canopy sublimation during simulation (m)
%     %     sumprec     1  0 summed precipitation during simulation (m)
%     %     sumsprec    1  0 summed snow precipitation during simulation (m)
%     %     sumunload   1  0 summed canopy unloading during the simulation (m)
%     %     sumroff     1  0 summed runoff during the simulation (m)
%     %     sumswemelt  1  0 summed snow-water-equivalent melt (m)
%     %     sumsublim   1  0 summed static-surface sublimation (m)
%     %     wbal        1  0 summed water balance error during the simulation (m)
%     if save_snowpack == 1
%         temp                =   squeeze(read_grads('snowpack.ctl','snowd'));
%         snowpack.snowd      = 	squeeze(temp(row,col,timestep)); clear temp
%         temp                =   squeeze(read_grads('snowpack.ctl','rosnow'));
%         snowpack.rosnow 	= 	squeeze(temp(row,col,timestep)); clear temp
%         temp                = 	squeeze(read_grads('snowpack.ctl','swed'));
%         snowpack.swed       =   squeeze(temp(row,col,timestep)); clear temp
%         temp                = 	squeeze(read_grads('snowpack.ctl','runoff'));
%         snowpack.runoff     =   squeeze(temp(row,col,timestep)); clear temp
%         temp                = 	squeeze(read_grads('snowpack.ctl','rain'));
%         snowpack.rain       =   squeeze(temp(row,col,timestep)); clear temp
%         temp                =   squeeze(read_grads('snowpack.ctl','sprec'));
%         snowpack.sprec      =   squeeze(temp(row,col,timestep)); clear temp
%         temp                = 	squeeze(read_grads('snowpack.ctl','sumqcs'));
%         snowpack.sumqcs 	=   squeeze(temp(row,col,timestep)); clear temp
%         temp                = 	squeeze(read_grads('snowpack.ctl','qcs'));
%         snowpack.qcs        =   squeeze(temp(row,col,timestep)); clear temp
%         temp                = 	squeeze(read_grads('snowpack.ctl','canopy'));
%         snowpack.canopy 	=   squeeze(temp(row,col,timestep)); clear temp
%         temp                = 	squeeze(read_grads('snowpack.ctl','sumqcs'));
%         snowpack.sumqcs 	=   squeeze(temp(row,col,timestep)); clear temp
%         temp                =   squeeze(read_grads('snowpack.ctl','sumprec'));
%         snowpack.sumprec 	= 	squeeze(temp(row,col,timestep)); clear temp
%         temp                =   squeeze(read_grads('snowpack.ctl','sumsprec'));
%         snowpack.sumsprec 	= 	squeeze(temp(row,col,timestep)); clear temp  
%         temp                = 	squeeze(read_grads('snowpack.ctl','sumunload'));
%         snowpack.sumunload  =   squeeze(temp(row,col,timestep)); clear temp
%         temp                = 	squeeze(read_grads('snowpack.ctl','sumroff'));
%         snowpack.sumrunoff  =   squeeze(temp(row,col,timestep)); clear temp
%         temp                = 	squeeze(read_grads('snowpack.ctl','sumswemelt'));
%         snowpack.sumswemelt =   squeeze(temp(row,col,timestep)); clear temp
%         temp                = 	squeeze(read_grads('snowpack.ctl','sumsublim'));
%         snowpack.sumsublim  =   squeeze(temp(row,col,timestep)); clear temp
%         temp                = 	squeeze(read_grads('snowpack.ctl','wbal'));
%         snowpack.error      =   squeeze(temp(row,col,timestep)); clear temp
%     end
%     
%     if delete_snowpack == 1
%         delete('snowpack.gdat');
%     end
% 
%     % ALBEDO
% 
%     if save_albedo == 1;
%         fname1 = 'albedo.ctl';
%         copyfile([ctlfold fname1],datafold);
%         albedo = squeeze(read_grads(fname1,'albedo'));
%         albedo = roundn(albedo,round_albedo);
%         albedo = rot90_3D(albedo,3,1);
%         data.albedo = squeeze(albedo(row,col,timestep));
%         clear albedo
%     end
% 
%     if delete_albedo == 1;
%         delete('albedo.gdat');
%     end
% 
%     % EBAL
% 
%     if save_ebal == 1;
%         fname2 = 'ebal.ctl';
%         copyfile([ctlfold fname2],datafold);
%         ebal = squeeze(read_grads(fname2,'ebal'));
%         ebal = roundn(ebal,round_ebal);
%         data.ebal = squeeze(ebal(row,col,timestep));
%         clear ebal
%     end
% 
%     if delete_ebal == 1;
%         delete('ebal.gdat');
%     end
%     
%     % GLACIER MELT
% 
%     if save_glacmelt == 1;
%         fname3 = 'glacmelt.ctl';
%         copyfile([ctlfold fname3],datafold);
%         glacmelt = squeeze(read_grads(fname3,'glacmelt'));
%         glacmelt = roundn(glacmelt,round_glacmelt);
%         data.glacmelt = squeeze(glacmelt(row,col,timestep));
%         clear glacmelt
%     end
% 
%     if delete_glacmelt == 1;
%         delete('glacmelt.gdat');
%     end
%     
%     % PREC
% 
%     if save_prec == 1;
%         fname4 = 'prec.ctl';
%         copyfile([ctlfold fname4],datafold);
%         prec = squeeze(read_grads(fname4,'prec'));
%         prec = roundn(prec,round_prec);
%         data.prec = squeeze(prec(row,col,timestep));
%         clear prec
%     end
% 
%     if delete_prec == 1;
%         delete('prec.gdat');
%     end
%     
%     % CONDUCTIVE HEAT FLUX
% 
%     if save_qc == 1;
%         fname5 = 'qc.ctl';
%         copyfile([ctlfold fname5],datafold);
%         qc = squeeze(read_grads(fname5,'qc'));
%         qc = roundn(qc,round_qc);
%         data.qc = squeeze(qc(row,col,timestep));
%         clear qc
%     end
% 
%     if delete_qc == 1;
%         delete('qc.gdat');
%     end
%     
%     % LATENT HEAT FLUX
% 
%     if save_qe == 1;
%         fname6 = 'qe.ctl';
%         copyfile([ctlfold fname6],datafold);
%         qe = squeeze(read_grads(fname6,'qe'));
%         qe = roundn(qe,round_qe);
%         data.qe = squeeze(qe(row,col,timestep));
%         clear qe
%     end
% 
%     if delete_qe == 1;
%         delete('qe.gdat');
%     end
%     
%     % REFREEZING HEAT FLUX
% 
%     if save_qf == 1;
%         fname7 = 'qf.ctl';
%         copyfile([ctlfold fname7],datafold);
%         qf = squeeze(read_grads(fname7,'qf'));
%         qf = roundn(qf,round_qf);
%         data.qf = squeeze(qf(row,col,timestep));
%         clear qf
%     end
% 
%     if delete_qf == 1;
%         delete('qf.gdat');
%     end
%     
%     % SENSIBLE HEAT FLUX
% 
%     if save_qh == 1;
%         fname8 = 'qh.ctl';
%         copyfile([ctlfold fname8],datafold);
%         qh = squeeze(read_grads(fname8,'qh'));
%         qh = roundn(qh,round_qh);
%         data.qh = squeeze(qh(row,col,timestep));
%         clear qh
%     end
% 
%     if delete_qh == 1;
%         delete('qh.gdat');
%     end
%     
%     % OUTGOING LONGWVE
% 
%     if save_qle == 1;
%         fname9 = 'qle.ctl';
%         copyfile([ctlfold fname9],datafold);
%         qle = squeeze(read_grads(fname9,'qle'));
%         qle = roundn(qle,round_qle);
%         data.qle = squeeze(qle(row,col,timestep));
%         clear qle
%     end
% 
%     if delete_qle == 1;
%         delete('qle.gdat');
%     end
%     
%     % INCOMING LONGWVE
% 
%     if save_qli == 1;
%         fname10 = 'qli.ctl';
%         copyfile([ctlfold fname10],datafold);
%         qli = squeeze(read_grads(fname10,'qli'));
%         qli = roundn(qli,round_qli);
%         data.qli = squeeze(qli(row,col,timestep));
%         clear qli
%     end
% 
%     if delete_qli == 1;
%         delete('qli.gdat');
%     end
%     
%     % MELT ENERGY
% 
%     if save_qm == 1;
%         fname11 = 'qm.ctl';
%         copyfile([ctlfold fname11],datafold);
%         qm = squeeze(read_grads(fname11,'qm'));
%         qm = roundn(qm,round_qm);
%         data.qm = squeeze(qm(row,col,timestep));
%         clear qm
%     end
% 
%     if delete_qm == 1;
%         delete('qm.gdat');
%     end
%     
%     % QSI
% 
%     if save_qsi == 1;
%         fname12 = 'qsi.ctl';
%         copyfile([ctlfold fname12],datafold);
%         qsi = squeeze(read_grads(fname12,'qsi'));
%         qsi = roundn(qsi,round_qsi);
%         data.qsi = squeeze(qsi(row,col,timestep));
%         clear qsi
%     end
% 
%     if delete_qsi == 1;
%         delete('qsi.gdat');
%     end
%     
%     % RAIN
% 
%     if save_rain == 1;
%         fname13 = 'rain.ctl';
%         copyfile([ctlfold fname13],datafold);
%         rain = squeeze(read_grads(fname13,'rain'));
%         rain = roundn(rain,round_rain);
%         data.rain = squeeze(rain(row,col,timestep));
%         clear rain
%     end
% 
%     if delete_rain == 1;
%         delete('rain.gdat');
%     end
%     
%     % RH
% 
%     if save_rh == 1;
%         fname14 = 'rh.ctl';
%         copyfile([ctlfold fname14],datafold);
%         rh = squeeze(read_grads(fname14,'rh'));
%         rh = roundn(rh,round_rh);
%         data.rh = squeeze(rh(row,col,timestep));
%         clear rh
%     end
% 
%     if delete_rh == 1;
%         delete('rh.gdat');
%     end
%     
%     % NEW SNOW DENSITY
% 
%     if save_ro_nsnow == 1;
%         fname15 = 'ro_nsnow.ctl';
%         copyfile([ctlfold fname15],datafold);
%         ro_nsnow = squeeze(read_grads(fname15,'ro_nsnow'));
%         ro_nsnow = roundn(ro_nsnow,round_ro_nsnow);
%         data.ro_nsnow = squeeze(ro_nsnow(row,col,timestep));
%         clear ro_nsnow
%     end
% 
%     if delete_ro_nsnow == 1;
%         delete('ro_nsnow.gdat');
%     end
%     
%     % SNOW DENSITY
% 
%     if save_rosnow == 1;
%         fname16 = 'rosnow.ctl';
%         copyfile([ctlfold fname16],datafold);
%         rosnow = squeeze(read_grads(fname16,'rosnow'));
%         rosnow = roundn(rosnow,round_rosnow);
%         data.rosnow = squeeze(rosnow(row,col,timestep));
%         clear rosnow
%     end
% 
%     if delete_rosnow == 1;
%         delete('rosnow.gdat');
%     end
%     
%     % RUNOFF
% 
%     if save_runoff == 1;
%         fname17 = 'runoff.ctl';
%         copyfile([ctlfold fname17],datafold);
%         runoff = squeeze(read_grads(fname17,'runoff'));
%         runoff = roundn(runoff,round_runoff);
%         data.runoff = squeeze(runoff(row,col,timestep));
%         clear runoff
%     end
% 
%     if delete_runoff == 1;
%         delete('runoff.gdat');
%     end
%     
%     % SNOW DEPTH
% 
%     if save_snowd == 1;
%         fname18 = 'snowd.ctl';
%         copyfile([ctlfold fname18],datafold);
%         snowd = squeeze(read_grads(fname18,'snowd'));
%         snowd = roundn(snowd,round_snowd);
%         data.snowd = squeeze(snowd(row,col,timestep));
%         clear snowd
%     end
% 
%     if delete_snowd == 1;
%         delete('snowd.gdat');
%     end
%     
%     % SPREC
% 
%     if save_sprec == 1;
%         fname19 = 'sprec.ctl';
%         copyfile([ctlfold fname19],datafold);
%         sprec = squeeze(read_grads(fname19,'sprec'));
%         sprec = roundn(sprec,round_sprec);
%         data.sprec = squeeze(sprec(row,col,timestep));
%         clear sprec
%     end
% 
%     if delete_sprec == 1;
%         delete('sprec.gdat');
%     end
% 
%     % SWE
% 
%     if save_swe == 1;
%         fname20 = 'swe.ctl';
%         copyfile([ctlfold fname20],datafold);
%         swe = squeeze(read_grads(fname20,'swe'));
%         swe = roundn(swe,round_swe);
%         data.swe = squeeze(swe(row,col,timestep));
%         clear swe
%     end
% 
%     if delete_swe == 1;
%         delete('swe.gdat');
%     end
%     
%     % SWEMELT
% 
%     if save_swemelt == 1;
%         fname21 = 'swemelt.ctl';
%         copyfile([ctlfold fname21],datafold);
%         swemelt = squeeze(read_grads(fname21,'swemelt'));
%         swemelt = roundn(swemelt,round_swemelt);
%         data.swemelt = squeeze(swemelt(row,col,timestep));
%         clear swemelt
% 
%     end
% 
%     if delete_swemelt == 1;
%         delete('swemelt.gdat');
%     end
%     
%     % TAIR. NOTE - need to add exist checks to rest of vars
%         
%     if save_tair == 1;
% 
%         fname1 = 'tair.ctl';
%         fname2 = [datafold 'tair.gdat'];
% 
%         if  ~exist(fname2,'file')
%             disp('File for saving not present in datafold. Skipping.')
%             missing_file_flag = 1;
%         else
%             missing_file_flag = 0;
%         end                
% 
%         if missing_file_flag == 0; % file exists
%             copyfile([ctlfold fname1],datafold);
%             tair = squeeze(read_grads(fname1,'tair'));                
%             naninds = find(isnan(tair)); % check for bad data
%             if ~isempty(naninds);
%                 disp(['NaN data detected, ' fname2]);
%                 pause
%             end
% 
%             tair = tair - 273.15;
%             tair = roundn(tair,round_tair); % no bad data, continue
% 
%             data.tair = squeeze(tair(row,col,timestep));
%             clear tair
% 
%             if delete_tair == 1; % delete the .gdat file
%                 delete('tair.gdat');
%             end
%         end
%     end
% 
%     % TMAX
% 
%     if save_tmax == 1;
%         fname1 = 'tair.ctl';
%         copyfile([ctlfold fname1],datafold);
%         tmax = squeeze(read_grads(fname1,'tair'));
%         tmax = tmax - 273.15;
%         tmax = roundn(tmax,round_tmax);
%         data.tmax = squeeze(tmax(row,col,timestep));
%         clear tmax
%     end
% 
%     if delete_tmax == 1;
%         delete('tair.gdat');
%     end
% 
%     % TMIN
% 
%     if save_tmin == 1;         
%         fname1 = 'tair.ctl';
%         copyfile([ctlfold fname1],datafold);
%         tmin = squeeze(read_grads(fname1,'tair'));
%         tmin = tmin - 273.15;
%         tmin = roundn(tmin,round_tmin);
%         data.tmin = squeeze(tmin(row,col,timestep));
%         clear tmin
%     end
% 
%     if delete_tmin == 1;
%         delete('tair.gdat');
%     end
% 
%     % SURFACE SKIN TEMP
%     
%     if save_tsfc == 1;         
%         fname23 = 'tsfc.ctl';
%         copyfile([ctlfold fname23],datafold);
%         tsfc = squeeze(read_grads(fname23,'tsfc'));
%         tsfc = roundn(tsfc,round_tsfc);
%         data.tsfc = squeeze(tsfc(row,col,timestep));
%         clear tsfc
%     end
% 
%     if delete_tsfc == 1;
%         delete('tsfc.gdat');
%     end
%     
%     % SUMMED WATER BALANCE
% 
%     if save_wbal == 1;
%         fname24 = 'wbal.ctl';
%         copyfile([ctlfold fname24],datafold);
%         wbal = squeeze(read_grads(fname24,'wbal'));
%         wbal = roundn(wbal,round_wbal);
%         data.wbal = squeeze(wbal(row,col,timestep));
%         clear wbal
%     end
% 
%     if delete_wbal == 1;
%         delete('wbal.gdat');
%     end
% 
%     % WDIR
%     if save_wdir == 1;
% 
%         fname7  = 'wdir.ctl';
%         fname72 = [datafold 'wdir.gdat'];
%         if  ~exist(fname72,'file')
%             disp([fname72 ' file not present in datafold. Skipping.'])
%             missing_file_flag = 1;
%         else
%             missing_file_flag = 0;
%         end                
% 
%         if missing_file_flag == 0; % file exists
%             copyfile([ctlfold fname7],datafold);
%             wdir = squeeze(read_grads(fname7,'wdir'));                
%             naninds = find(isnan(wdir)); % check for bad data
%             if ~isempty(naninds);
%                 disp(['NaN data detected, ' fname2]);
%                 pause
%             end
% 
%             wdir = roundn(wdir,round_wdir);
%             data.wdir = squeeze(wdir(row,col,timestep));
%             clear wdir
% 
%             if delete_wdir == 1;
%                 delete('wdir.gdat');
%             end
%         end                    
%     end
% 
%     % WSPD
% 
%     % wspd
%    if save_wspd == 1;
% 
%         fname8  = 'wspd.ctl';
%         fname82 = [datafold 'wspd.gdat'];
%         if  ~exist(fname82,'file')
%             disp([fname82 ' file not present in datafold. Skipping.'])
%             missing_file_flag = 1;
%         else
%             missing_file_flag = 0;
%         end                
% 
%         if missing_file_flag == 0; % file exists
%             copyfile([ctlfold fname8],datafold);
%             wspd = squeeze(read_grads(fname8,'wspd'));                
%             naninds = find(isnan(wspd)); % check for bad data
%             if ~isempty(naninds);
%                 disp(['NaN data detected, ' fname2]);
%                 pause
%             end
% 
%             wspd = roundn(wspd,round_wspd);
%             data.wspd = squeeze(wspd(row,col,timestep));
%             clear wspd
% 
%             if delete_wspd == 1;
%                 delete('wspd.gdat');
%             end
%         end                    
%     end
%        
%     %% Save the data
%     if save_enbal == 1
%         save([savefold 'enbal'],'enbal');
%         clear enbal
%     end
%     if save_micromet == 1
%         save([savefold 'micromet'],'micromet');
%         clear micromet
%     end
%     if save_snowpack == 1
%         save([savefold 'snowpack'],'snowpack');
%         clear snowpack
%     end
%     if save_data == 1
%         save([savefold 'data'],'data');
%         clear data
%     end
% end
% 
%     
     