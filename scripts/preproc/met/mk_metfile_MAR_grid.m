clean

% NOTE: when I ran this on my personal mac I should have used xmask/ymask
% not xmodel/ymodel to get the runs on the same numbering as the original
% runs which was xmask/ymask. the mistake was in using list=dir() in THIS
% script to save the data, now that mistake is propagated into the new runs
% so i'll just have to keep using xmodel/ymodel

savedata    =  false;
sitename    =  'region';
yri         =  2008;
yrf         =  2018;
nyears      =  yrf-yri+1;
tinterp     =  true;
newdt       =  '15m';

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% set paths
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pathdata    =  setpath('mar3.11','data');
pathsave    =  setpath('mar3.11/matfiles/region/level2/grids/','data');
pathmodis   =  setpath('geus/albedo/gris/','data');
filelist    =  getlist(pathdata,'*.nc');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% read in the catchment to get the x,y coordinate 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% load data
load('projsipsn.mat','projsipsn')         % sea ice polar stereographic north
load('modis_ice_mask.mat','icemask');     % bare ice mask x,y coordinates
load('imbie_SW.mat','imbie');             % imbie sectors

% use the bare ice mask x,y coord's to build sector met files
Xmet     = icemask.xmodel;
Ymet     = icemask.ymodel;

% use the sw sector as an ROI mask
ROI      = imbie.poly;

% some variables need only be extracted once, do them first
fmar     = [pathdata 'MARv3.11-ERA5-15km-2018.nc'];
info     = ncparse(fmar);
LON      = ncread(fmar,'LON');
LAT      = ncread(fmar,'LAT');
VAR      = ncread(fmar,'RUH');
[Xmar,...
   Ymar] = projfwd(projsipsn,LAT,LON);


% figure; scatter(Xmar(:),Ymar(:)); hold on; scatter(Xmet,Ymet);

% similarly, do the modis conversion from lat/lon to x/y first
fmodis   = [pathmodis 'Greenland_Reflectivity_' num2str(yri) '_5km_C6.nc'];
LON      = double(ncread(fmodis,'lon'));
LAT      = double(ncread(fmodis,'lat'));
[Xmod, ...
   Ymod] = projfwd(projsipsn,LAT,LON);

% for testing:
n=nyears;

for n = 1:nyears
    
   yyyy    = yri+n-1;
   fmar    = [pathdata 'MARv3.11-ERA5-15km-' num2str(yyyy) '.nc'];
   fmodis  = [pathmodis 'Greenland_Reflectivity_' num2str(yyyy) '_5km_C6.nc'];
   
   % for the region gridded runs, use the simple naming convention
   fsavemet    = [pathsave 'met_' int2str(yyyy) '.mat'];
   
%    commented out for testing
%    if ~exist(savepath,'dir'); mkdir(savepath); end
      
   % the filename gets appended with the point number in the function
   fsavedata   =  [pathsave 'data_' int2str(yyyy) '.mat'];

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
   [met,data]  =  makeMarMetfileGrid(fmar,'Xmar',Xmar,'Ymar',Ymar,      ...
                  'ROI',ROI,'Xmet',Xmet,'Ymet',Ymet,'fmodis',fmodis,    ...
                  'Xmodis',Xmod,'Ymodis',Ymod);
               

   [met,data]  =  makeMarMetfileGrid(fmar,'Xmar',Xmar,'Ymar',Ymar,      ...
                  'ROI',ROI,'fmodis',fmodis,'Xmodis',Xmod,'Ymodis',Ymod);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if savedata == true
      save(fsavemet,'met','-v7.3')
      save(fsavedata,'data','-v7.3')
   end

   clear met data
   
end


% test the met file
% vars    =   met.Properties.VariableNames;
% units   =   met.Properties.VariableUnits;

    % % this will be useful for generating files for every grid cell        
    %         latstr  = sprintf('%.6f',met.Lat(1));
    %         lonstr  = strrep(sprintf('%.6f',met.Lon(1)),'-','');
    %         fsave   = [fsave  '_' latstr '_' lonstr '.mat'];
    %         save(fsave,'met');
    
% % PICK UP HERE = basically ready to build a met file, see loop below,
% % consider renaming vars
% 
% % using a polygon
% polyb   = polybuffer(poly,12000);
% xpolyb  = polyb.Vertices(:,1);
% ypolyb  = polyb.Vertices(:,2);
% 
% [xstart,xcount] =   ncrowcol(X,X,Y,xpolyb,ypolyb);
% [ystart,ycount] =   ncrowcol(Y,X,Y,xpolyb,ypolyb);
% [vstart,vcount]     =   ncrowcol(squeeze(VAR(:,:,:,1)),X,Y,xpolyb,ypolyb);
% 
% % add time start/end to srof,crof with a daily calendar the data is stored
% % as row,col,1:24,1:365 
% tdaily           = datetime(yyyy,1,1):days:datetime(yyyy,12,31);
% [si,ei]     = dateInds(datetime(yyyy,6,1),datetime(yyyy,9,1),tdaily);
% vstart        = [vstart si];          
% vcount        = [vcount (ei-si)];
% 
% 
% % this was in compare_racmo_versions
% 
% % tried this to use ncrowcol before realizing it's not needed with my point
% % interpolation approach
% % get the lat/lon indices to read in the data around the catchment
% % xpolymin    =   polymin.Vertices(:,1);
% % xpolymax    =   polymax.Vertices(:,1);
% % [slat,clat] =   ncrowcol(LAT,rot90(X,3),rot90(Y,3),xpolyb,ypolyb);
% % [slon,clon] =   ncrowcol(LON,rot90(X,3),rot90(Y,3),xpolyb,ypolyb);
% % [sswd,cswd] =   ncrowcol(swsd,rot90(X,3),rot90(Y,3),xpolyb,ypolyb);
% 
% % ... maybe I can just use the .nc files with the same workflow ...
% 
% 
% 
% 
% 
% 
% % here's a way it could be done, from test_RACMO_refreeze
% 
% vers    = 'subsurface';
% 
% % set paths, get files
% pathdata  = ['/Users/coop558/mydata/racmo2.3/' vers '/'];
% p.mask  = ['/Users/coop558/mydata/racmo2.3/' vers '/'];
% list    = dir(fullfile([pathdata '*.nc']));
% 
% load('projsipsn.mat');
% load('racmo_subsurf_time.mat');
% 
% fnames      =   {list.name};
% info        =   ncparse([pathdata list(1).name]);
% lon         =   readRacmo2p3([pathdata list(1).name],'lon');
% lat         =   readRacmo2p3([pathdata list(1).name],'lat');
% 
% % is runoff ~= meltin + snowmelt - refreeze?
% % get coordinates in the southwest
% latdif      =   abs(lat-65);
% londif      =   abs(lon+48);
% [r,c]       =   find(latdif==min(latdif) & londif==min(londif));
% t1          =   find(T==datetime(2016,7,1));
% t2          =   find(T==datetime(2016,9,1));
% tdaily           =   T(t1:t2);
% ameltin     =   squeeze(meltin(r,c,t1:t2));
% asnowmelt   =   squeeze(snowmelt(r,c,t1:t2));
% arunoff     =   squeeze(runoff(r,c,t1:t2));
% afreeze     =   squeeze(freeze(r,c,t1:t2));
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % BELOW here is from the build metfiles script, above is testing to figure
% % out how to do it directly from teh .nc file using ncrowcol
% 
% 
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %% paths
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% pathdata      =   '/Users/coop558/mydata/mar3.11/matfiles/region/level1/';
% pathsave      =   ['/Users/coop558/mydata/mar3.11/matfiles/region/level2/'...
%                     albedo '/'];
% pathsupp     =   ['/Users/coop558/mydata/DATA/greenland/geus/albedo/'    ...
%                     'region/matfiles/' proj '/native/'];
% list        =   dir(fullfile([pathdata '*.mat']));
% 
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                            
% %% load the modis ice mask, extract x,y locations
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% load([pathsupp 'mask/modis_ice_mask']);
% 
% mask        =   icemask.mask;
% elev        =   icemask.elev;
% x           =   icemask.x;
% y           =   icemask.y;
% [nr,nc]     =   size(mask);
% maskrs      =   logical(reshape(mask,nr*nc,1));
% 
% % reshape the elevation and x/y coord's and apply the mask
% elev        =   reshape(elev,nr*nc,1);
% x           =   reshape(x,nr*nc,1);
% y           =   reshape(y,nr*nc,1);
% x           =   x(maskrs);
% y           =   y(maskrs);
% elev        =   elev(maskrs);
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %% use this to figure out the indices for rio behar (it's 1045)
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % load('catchment_2016.mat');
% % figure; 
% % scatter(x,y,40,elev,'filled'); colorbar; hold on;
% % plot(catchment.bestguess_x,catchment.bestguess_y);
% % idx = find(x==-1.725e5 & y==-2.5125e6);
% % [r,c] = ind2sub(size(icemask.mask),idx);
% % scatter(x(idx),y(idx),400,'filled')
% % plot(catchment.bestguess_x,catchment.bestguess_y);
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %% loop through all years and build metfiles for all the points
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% for n = 1:nyears
%     
%     yyyy              =   yri+n-1;
%     fname           =   [pathdata 'mar_' int2str(yyyy)];
%     load(fname)
%     
% %--------------------------------------------------------------------------
%     % first read in modis albedo if requested and convert to hourly
% %--------------------------------------------------------------------------    
% 
%   % build a daily calendar 
%     startdate           =   [yyyy,1,1,12,0,0];
%     enddate             =   [yyyy,12,31,12,0,0];
%     t_daily             =   (datenum(startdate):1:datenum(enddate))';
%     
%   % build a 1-hourly calendar on to which the daily data will be resampled
%     startdate           =   [yyyy,1,1,0,0,0];
%     enddate             =   [yyyy,12,31,23,0,0];
%     t_1hrly             =   (datenum(startdate):1/24:datenum(enddate))';
%     nhrs                =   length(t_1hrly);
% 
%   % read in the modis data and interpolate it to 1-hourly
%     if strcmp(albedo,'modis')
%         load([pathsupp 'modis_albedo_' int2str(yyyy)]);
%         [nr,nc,nd]  =   size(modis.A);
%         alb         =   reshape(modis.A,nr*nc,nd);
%         alb         =   (alb(maskrs,:))';
%         if ~isleap(yyyy)
%             alb     =   alb(1:365,:);
%         end
%         alb         =   interp1(t_daily,alb,t_1hrly,'linear');
%     else
%         alb         =   (mar.albedo(maskrs,:))';
%     end
%     
% %--------------------------------------------------------------------------
%     % deal with the remaining variables
% %--------------------------------------------------------------------------   
% 
%   % pull out the relevant forcing data and apply the mask 
%     [Y,M,D,H,~,~]   =   datevec(mar.dates);
%     swd             =   (mar.swd(maskrs,:))';
%     lwd             =   (mar.lwd(maskrs,:))';
%     tair            =   (mar.tair(maskrs,:))';
%     relh            =   (mar.rh(maskrs,:))';
%     wspd            =   (mar.wspd(maskrs,:))';
%     wdir            =   (mar.wdir(maskrs,:))';
%     ppt             =   (mar.prec(maskrs,:))';
%     press           =   (mar.psfc(maskrs,:))';
%     [nhrs,npts]     =   size(swd);
% 
%   % make the met files
%     for n = 1:npts
% 
%         id_n        =   n.*ones(nhrs,1);
%         x_n         =   x(n).*ones(nhrs,1);
%         y_n         =   y(n).*ones(nhrs,1);
%         elev_n      =   roundn(elev(n),0).*ones(nhrs,1);
%         tair_n      =   linearfill(roundn(tair(:,n),-2));
%         relh_n      =   linearfill(roundn(relh(:,n),0));
%         wspd_n      =   linearfill(roundn(wspd(:,n),-1));
%         wdir_n      =   linearfill(roundn(wdir(:,n),0));
%         ppt_n       =   linearfill(roundn(ppt(:,n),-3));
%         swd_n       =   linearfill(roundn(swd(:,n),-1));
%         lwd_n       =   linearfill(roundn(lwd(:,n),-1));
%         alb_n       =   linearfill(roundn(alb(:,n),-3));
%         press_n     =   linearfill(roundn(press(:,n),-3));
%         
%         met         =   [Y,M,D,H,id_n,x_n,y_n,elev_n,tair_n,relh_n,     ...
%                         wspd_n,wdir_n,ppt_n,swd_n,lwd_n,alb_n,press_n];
% 
%         if save_data == 1
%             psave   =   [pathsave int2str(yyyy)];
%             if ~exist(psave,'dir')
%                 mkdir(psave)
%             end
%             fsave   =   [psave '/met_' int2str(id_n(1)) '.mat'];
%             save(fsave,'met')
%         end
%         clear met
%         
%     end
%     
% end
% 
% % make figures to confirm things make sense
% % if plot_figs == 1
% % figure
% % scatter(x,y,20,swd(4000,:),'filled')
% % title('SW Down'); colorbar
% % 
% % figure
% % scatter(x,y,20,lwd(4000,:),'filled')
% % title('LW Down'); colorbar
% % 
% % figure
% % scatter(x,y,20,elev,'filled')
% % title('Elevation'); colorbar
% % end
% 
% 
% 
% % % I shoudl be able to get RH from this:
% % Tf  =   273.16;
% % A   =   6.1115 * 100.0;
% % B   =   22.452;
% % C   =   272.55;
% % % Compute the saturated water vapor pressure at the surface.
% % es0 =   A.*exp((B.*(mar.TTH - Tf))./(C + (mar.TTH - Tf))); 
% % es  =   mar.QQH.*461.*mar.TTH;
% % RH  =   es./es0;
