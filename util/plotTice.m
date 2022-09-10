
% this all came from v8/hills/c_icemodel

function h = plotTice(opts,ice2)

   site     = opts.sitename;
   yyyy     = opts.yyyy;
   deltaz   = opts.dz_thermal;
   yr       = str2double(yyyy);

   % paths to the ice temperature and ice density data
   patheval = setpath('GREENLAND/runoff/data/icemodel/eval/');

   switch yyyy
      case '2014'
         load([patheval 'hills/hills_Tice_' yyyy],'T14'); T = T14; clear T14;
      case '2015'
         load([patheval 'hills/hills_Tice_' yyyy],'T15'); T = T15; clear T15;
      case '2016'
         load([patheval 'hills/hills_Tice_' yyyy],'T16'); T = T16; clear T16;
   end

   load([patheval '2016/density_2016.mat'],'density'); rho16=density;
   load([patheval '2016/density_2018.mat'],'density'); rho18=density;

   % Calculate ice temperature
   Z        =  ice2.Z;
   Time     =  ice2.Time;
   t1       =  datetime(yr,7,1,0,0,0,'TimeZone','UTC');
   t2       =  datetime(yr,7,31,23,0,0,'TimeZone','UTC');
   idx      =  isbetween(Time,t1,t2);
   Tice     =  ice2.Tice(:,idx);
   ndays    =  size(Tice,2)/24;
   nz       =  size(Tice,1);
   Tice     =  reshape(Tice,nz,24,ndays);
   TiceD    =  mean(squeeze(mean(Tice(:,12:14,:),2)),2);
   TiceN    =  mean(squeeze(mean(Tice(:,2:4,:),2)),2);

   % Ben Hill's temperature data
   di    =  1;
   if isfield(T,'depth')
      T.Depth = T.depth;
   end
     
   de =  find(T.Depth(:,1) == deltaz*nz);
   
   if isempty(de)
       de = size(T.Y,1);
   end
   val.iceT       =  T.T(di:de,idx);
   val.iceT       =  reshape(val.iceT,size(val.iceT,1),24,ndays);
   val.iceTday    =  mean(squeeze(mean(val.iceT(:,12:14,:),2)),2);
   val.iceTnight  =  mean(squeeze(mean(val.iceT(:,2:4,:),2)),2);
   val.z          =  (T.Depth(di:de))';
   zi2            =  find(round(val.z,2) == 2.0);
   zi             =  find(round(Z,1)==2,1,'last');
   
   % ltext = {sprintf('\nMeasured\n(Hills et al. 2018)'),'IceModel'};
   ltext = {'Measured','(Hills et al. 2018)','IceModel'};
   dtext = {'Daytime Summer Ice Temperature ($^\circ$C)','Depth (m)'};
   ntext = {'Nighttime Summer Ice Temperature ($^\circ$C)','Depth (m)'};


   % ICE TEMPERATURE DAY TIME
   macfig;
   subplot(2,2,1);
   plot(val.iceTday(1:zi2),val.z(1:zi2),'o'); hold on; 
   plot(1,1,'Color','none'); plot(TiceD(1:zi),ice2.Z(1:zi)); 
   set(gca,'YDir','reverse'); 
   legend(ltext,'Location','best'); xylabel(dtext{:});
   formatPlotMarkers('markersize',6)

   % ICE TEMPERATURE NIGHT TIME
   subplot(2,2,2);
   plot(val.iceTnight(1:zi2),val.z(1:zi2),'o'); hold on; 
   plot(1,1,'Color','none'); plot(TiceN(1:zi),ice2.Z(1:zi));
   set(gca,'YDir','reverse');
   legend(ltext,'Location','best'); xylabel(ntext{:});
   formatPlotMarkers('markersize',6)

   % DEEP ICE TEMPERATRUE
   subplot(2,2,3);
   plot(val.iceTday,val.z,'o'); hold on; 
   plot(1,1,'Color','none'); plot(TiceD,ice2.Z); 
   set(gca,'YDir','reverse'); 
   legend(ltext,'Location','best'); xylabel(dtext{:});
   formatPlotMarkers('markersize',6)

   % FIGURE 20 = ICE TEMPERATURE NIGHT TIME
   subplot(2,2,4);
   plot(val.iceTnight,val.z,'o'); hold on; 
   plot(1,1,'Color','none'); plot(TiceN,ice2.Z); 
   set(gca,'YDir','reverse');
   legend(ltext,'Location','best'); xylabel(ntext{:});
   formatPlotMarkers('markersize',6)



   % MEAN ANNUAL ICE TEMPERATURE
   Tobs  = mean(T.T(di:de,:),2);
   Zobs  = (T.Depth(di:de))';
   Tmod  = mean(ice2.Tice,2);
   Zmod  = ice2.Z;

   macfig('size','large');
   plot(Tobs,Zobs,'o'); hold on; plot(Tobs(1),Zobs(1),'Color','none');
   plot(Tmod,Zmod); 
   set(gca,'YDir','reverse'); 
   legend(ltext,'Location','best'); xylabel(dtext{:});
   formatPlotMarkers('markersize',6)
   title('mean annual ice temperature');


   %%%%%%%%%

%    % Ground heat flux
%    GHF             =   ice2.Qg(:,si:ei-1);
%    ndays           =   size(GHF,2)/24;
%    GHF             =   reshape(GHF,JJ,24,ndays);
%    GHFday          =   mean(squeeze(GHF(:,13,:)),2);
%    GHFnight        =   mean(squeeze(GHF(:,3,:)),2);

   %%%%%%%%%

% % these were all in a section 'if skinmodel false'
%    ice2.ro_ice    = ice2.f_ice.*917;
%    ice2.rel_sat   = ice2.f_liq ./ (1-ice2.f_ice);
% 
%    % MESH ICE TEMPERATURE 
%    plotice2(ice2,'Tice','hoursofday',1:24);
% 
%    % MESH ICE DENSITY
%    plotice2(ice2,'ro_ice','hoursofday',1:24);
% 
%    % MESH WATER FRAC
%    plotice2(ice2,'f_liq','hoursofday',1:24);
% 
%    % MESH REL SAT
%    plotice2(ice2,'f_liq','hoursofday',1:24);




%    % SUBPLOT version of Ground Heat Flux
% 
%    f18 = figure('Units','inches','Position',[5 5 8 8]);
%    subplot(1,2,1)
%    plot(10.*GHFday(1:i1),ice2.Z(1:i1));
%    set(gca,'YDir','reverse'); 
%    ylabel('Depth (m)');
%    xlabel('Summer GHF (W m^{-2})');
%    legend('13:00 UTC');
%    box off
% 
%    subplot(1,2,2)
%    plot(10.*GHFnight(2:i1),ice2.Z(2:i1)); 
%    set(gca,'YDir','reverse'); 
%    ylabel('Depth (m)');
%    xlabel('Summer GHF (W m^{-2})');
%    legend('03:00 UTC');
%    box off


%    % FIGURE 21 = GHF DAY vs NIGHT
%    f21 = figure('Units','inches','Position',[5 5 5 8]);
%    plot(GHFday(1:zi),ice2.Z(1:zi)); hold on;
%    plot(GHFnight(1:zi),ice2.Z(1:zi)); 
%    set(gca,'YDir','reverse'); 
%    ylabel('Depth (m)');
%    xlabel('Summer GHF (W m^{-2})');
%    legend('13:00 UTC','03:00 UTC','location','best');
%    box off
% 
% % FIGURE 24 = MESH GROUND HEAT FLUX
% si                  =   100; % doy 200
% ei                  =   300; % end of year
% Z                   =   ice2.Qg;
% Zrs                 =   reshape(Z,JJ,24,size(Z,2)/24);
% Z                   =   squeeze(mean(Zrs,2));
% f8                  =   figure;
% h                   =   pcolor(X,Y,Z(1:zi,si:ei)); ax = gca; datetick('x');
% h.EdgeColor         =   'none';
% ax.YDir             =   'reverse';
% ax.YLabel.String    =   'Depth (m)';
% ax.YLim             =   [0 2];
% c                   =   colorbar;
% axis tight
% title('GROUND HEAT FLUX (W m^{-2})')
% 
% % FIGURE 25 = GROUND HEAT FLUX AT 13:00
% Z                   =   ice2.Qg;
% Zrs                 =   reshape(Z,JJ,24,size(Z,2)/24);
% Z1                  =   squeeze(Zrs(:,13,:));
% Z2                  =   squeeze(Zrs(:,3,:));
% 
% f25                 =   figure;
% h                   =   pcolor(X,Y,Z1(1:zi,si:ei)); ax = gca; datetick('x');
% h.EdgeColor         =   'none';
% ax.YDir             =   'reverse';
% ax.YLabel.String    =   'Depth (m)';
% ax.YLim             =   [0 2];
% c                   =   colorbar;
% axis tight
% title('GROUND HEAT FLUX (W m^{-2}) 13:00 UTC')
% 
% % FIGURE 26 = GROUND HEAT FLUX AT 03:00
% f26                 =   figure;
% h                   =   pcolor(X,Y,Z2(1:zi,si:ei)); ax = gca; datetick('x');
% h.EdgeColor         =   'none';
% ax.YDir             =   'reverse';
% ax.YLabel.String    =   'Depth (m)';
% ax.YLim             =   [0 2];
% c                   =   colorbar;
% axis tight
% title('GROUND HEAT FLUX (W m^{-2}) 03:00 UTC')



