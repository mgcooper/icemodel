clean

% given a point, this builds a mar metfile. this limits the current
% function to regions that are well-approximated by a single grid cell.

% jul 2022 update, removed the saving in two places, reconcicled the naming
% with the filenames in icemodel/input/met, added modis albedo by default
% i.e. no more making two separate files just for the modis albedo

savedata    =  true;
modis       =  true;
sitename    =  'hills';
yri         =  2014;
yrf         =  2016;
nyears      =  yrf-yri+1;
tinterp     =  true;
newdt       =  '15m';

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% set paths
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% starting with hills data, I am only saving in the icemodel input folder
pathdata    =  '/Users/coop558/mydata/mar3.11/RUH2/';
pathsave    =  setpath('GREENLAND/icemodel/input/met/');
pathsupp    =  '/Users/coop558/mydata/geus/albedo/raw/';
list        =  getlist(pathdata,'*.nc');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% read in the catchment to get the x,y coordinate 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load('projsipsn.mat')

if strcmp(sitename,'hills')
   load('hills_Tice.mat');
   Xmet     = temperature.T_14.x;
   Ymet     = temperature.T_14.y;
else
   load(['Q_' sitename '.mat'])
   if isfield(sf,'med')
       Xmet    = nanmedian(sf.med.spsn.x);
       Ymet    = nanmedian(sf.med.spsn.y);
   elseif isfield(sf,'max')
       Xmet    = nanmedian(sf.max.spsn.x);
       Ymet    = nanmedian(sf.max.spsn.y);
   end
end


% some variables need only be extracted once, do them first
fmar        = [pathdata 'MARv3.11-ERA5-15km-' num2str(yri) '.nc'];
info        = ncparse(fmar);
LON         = ncread(fmar,'LON');
LAT         = ncread(fmar,'LAT');
VAR         = ncread(fmar,'RUH');
[Xmar,Ymar] = projfwd(projsipsn,LAT,LON);


for n = 1:nyears
    
   yyyy    = yri+n-1;
   fmar    = [pathdata 'MARv3.11-ERA5-15km-' num2str(yyyy) '.nc'];
   fmodis  = [pathsupp 'Greenland_Reflectivity_' num2str(yyyy) '_5km_C6.nc'];
   fsave   = ['met_' sitename '_MAR_' num2str(yyyy)];

% % we use 2014 for each year but keep this in case we use the 15/16 or 11
% data (2014 is the install date, the data is useable for 2015/16)
%    if strcmp(sitename,'hills')
%       if yyyy == 2014
%          Xmet  = temperature.T_14.x;
%          Ymet  = temperature.T_14.y;
%       elseif yyyy == 2015
%          Xmet  = (temperature.T_15a.x + temperature.T_15b.x + ...
%                      temperature.T_15c.x)/3;
%          Ymet  = (temperature.T_15a.y + temperature.T_15b.y + ...
%                      temperature.T_15c.y)/3;
%       elseif yyyy == 2016
%          Xmet  = temperature.T_16.x;
%          Ymet  = temperature.T_16.y;
%       end
%    end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
   met = makeMarMetfile(fmar,Xmar,Ymar,Xmet,Ymet,yyyy, ...
                           'modis',modis,'fmodis',fmodis);
        
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if savedata == true
      save([pathsave fsave '_1hr.mat'],'met');

      % interpolate, re-check out-of-bound values, and save
      if tinterp == true
         met      = interpMet(met,newdt);
         met.date = datenum(met.Time);
         met      = metchecks(met,false); % false = don't plot

         save([pathsave fsave '_' newdt '.mat'],'met');
      end
   end
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
