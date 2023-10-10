clean

% given a point, this builds a mar metfile. need to add option to build
% files for a given number of surrounding points, and/or all points within
% a polygon, and/or interpolation of points. this limits the current
% function to regions that are well-approximated by a single grid cell.

% note: mar RUH has problematic pixesl

% It would be helpful to have snow depth and skin temperature

save_data   =   true;
plot_check  =   true;
albedo      =   {'mar','modis'};      % 'modis' or 'mar'
basin       =   'upper_basin';
yri         =   2016;
yrf         =   2016;
nyears      =   yrf-yri+1;

testbuffer  =   false;                  % set false once buffer is correct 
res         =   50;                    % resampling resolution
buffer      =   14000;                   % this requires guess and check. 
method      =   'natural';              % resampling method

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% set paths
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p.mar       =   '/Users/coop558/mydata/mar3.11/RUH2/';
p.save      =   ['/Users/coop558/mydata/mar3.11/matfiles/' basin '/met/'];
p.copy      =   setpath('GREENLAND/icemodel/model/input/met/');
p.modis     =   '/Users/coop558/mydata/geus/albedo/raw/';
list        =   getlist(p.mar,'*.nc');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% read in the catchment to get the x,y coordinate 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% this can be used to get a single RCM point for a catchment, or if i add
% functinality to capture a grid of points, it can be used to interpolate
% them onto the catchment 

load(['Q_' basin '.mat'])
if isfield(sf,'med')
    Xmet    = nanmedian(sf.med.spsn.x);
    Ymet    = nanmedian(sf.med.spsn.y);
elseif isfield(sf,'max')
    Xmet    = nanmedian(sf.max.spsn.x);
    Ymet    = nanmedian(sf.max.spsn.y);
end

load('projsipsn.mat')
        
% pull out the lat lon and orient it correcly
fname       = [p.mar 'MARv3.11-ERA5-15km-' num2str(yri) '.nc'];
info        = ncparse(fname);
LON         = double(flipud(permute(ncread(fname,'LON'),[2 1])));
LAT         = double(flipud(permute(ncread(fname,'LAT'),[2 1])));
[X,Y]       = projfwd(projsipsn,LAT,LON);
T           = datetime(yri,1,1,0,0,0):hours(1):datetime(yri,12,31,23,0,0);

% this works for hourly variable (4-d) or daily (3-d)
% mass is in mmWE/hr, convert to m/day for an intuitive magnitude
VAR         = (24/1000).*double(ncread(fname,'RUH'));
[n,m,h,d]   = size(VAR);            % nrows x mcols x hours x days
VAR         = reshape(VAR,n,m,h*d);
VAR         = flipud(permute(VAR,[2 1 3]));

% mar has no-data values that are either 999 or 1e36
VAR         = setnan(VAR,[],VAR>998);
is          = find(T==datetime(yri,5,1,0,0,0));
ie          = find(T==datetime(yri,9,30,0,0,0));
Z           = VAR(:,:,is:ie);
max(Z(:))

worldmap('Greenland'); 
scatterm(LAT(:),LON(:),12,tocol(mean(Z,3)),'filled');
colorbar; title('RUH')

% for testing
n = 1; m = 1;
fmodis = [p.modis 'Greenland_Reflectivity_' num2str(yri) '_5km_C6.nc'];

% options 
opts.bufferDistance =   buffer;
opts.testBuffer     =   testbuffer;  
opts.siteName       =   basin;
opts.newResolution  =   res;
opts.method         =   method;
opts.saveFile       =   save_data;
opts.useModis       =   true;
opts.modisFilename  =   fmodis;

poly    =   sf.med.spsn.poly;
X       =   X(:);
Y       =   Y(:);
V       =   tocol(mean(Z,3));

data    =   saveMarData(fname,X,Y,poly,opts,V);


% if plot_check
%     hplot   = plot_grid(X,Y,LAT,LON,VAR);
% end
% 
% for n = 1:nyears
%     
%     yyyy    = yri+n-1;
%     fname    = [p.mar 'MARv3.11-ERA5-15km-' num2str(yyyy) '.nc'];
%     fmodis  = [p.modis 'Greenland_Reflectivity_' num2str(yyyy) '_5km_C6.nc'];
%     fsave   = ['test_' basin '_' num2str(yyyy)];
%     
%     for m = 1:numel(albedo)
%         
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%         
%       met = mk_mar_metfile(fname,X,Y,Xmet,Ymet,yyyy,albedo{m},fmodis);
%         
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
%         if save_data == true
%             save([p.save fsave '_' albedo{m} '_1hr'],'met');
%             save([p.copy fsave '_' albedo{m} '_1hr'],'met');
% 
%             % interpolate, re-check out-of-bound values, and save
%             if tinterp == true
%                 met      = interpmet(met,newdt);
%                 met.date = datenum(met.Time);
%                 met      = metchecks(met,false); % false = don't plot
%                 
%                 save([p.save fsave '_' albedo{m} '_' newdt],'met');
%                 save([p.copy fsave '_' albedo{m} '_' newdt],'met');
%             end
%         end
%     end
% end
% 
% 
% function h = plot_grid(X,Y,LAT,LON,VAR,plot_check)
%     
% VAR   = mean(VAR,3,'omitnan');
% 
% % with these, geoshow(LAT,LON) won't work, b/c it needs lat/lon grids
% % VAR = VAR(:);
% % X   = X(:);
% % Y   = Y(:);
% % LAT = LAT(:);
% % LON = LON(:);
% 
% % plot the average runoff using the x/y and lat/lon    
%     macfig;
%     
%     subplot(2,2,1);
%     scatter(X(:),Y(:),20,VAR(:),'filled'); colorbar; 
%     title('scatter(X,Y)');
%     
%     subplot(2,2,2);
%     scatter(LON(:),LAT(:),20,VAR(:),'filled'); colorbar; 
%     title('scatter(LON,LAT)');
%     
%     subplot(2,2,3);
%     worldmap('Greenland'); scatterm(LAT(:),LON(:),12,VAR(:),'filled'); colorbar
%     title('scatterm(LAT,LON)');
%     
%     subplot(2,2,4);
%     geoshow(LAT,LON,VAR,'DisplayType','texturemap');
%     title('geoshow(LAT,LON)');
% 
% h.gcf   = gcf;
% h.gca   = gca;
% 
% end
% 
% % test the met file
% % vars    =   met.Properties.VariableNames;
% % units   =   met.Properties.VariableUnits;
% 
%     % % this will be useful for generating files for every grid cell        
%     %         latstr  = sprintf('%.6f',met.Lat(1));
%     %         lonstr  = strrep(sprintf('%.6f',met.Lon(1)),'-','');
%     %         fsave   = [fsave  '_' latstr '_' lonstr '.mat'];
%     %         save(fsave,'met');
%     
% % % PICK UP HERE = basically ready to build a met file, see loop below,
% % % consider renaming vars
% % 
% % % using a polygon
% % polyb   = polybuffer(poly,12000);
% % xpolyb  = polyb.Vertices(:,1);
% % ypolyb  = polyb.Vertices(:,2);
% % 
% % [xstart,xcount] =   ncrowcol(X,X,Y,xpolyb,ypolyb);
% % [ystart,ycount] =   ncrowcol(Y,X,Y,xpolyb,ypolyb);
% % [vstart,vcount]     =   ncrowcol(squeeze(VAR(:,:,:,1)),X,Y,xpolyb,ypolyb);
% % 
% % % add time start/end to srof,crof with a daily calendar the data is stored
% % % as row,col,1:24,1:365 
% % tdaily           = datetime(yyyy,1,1):days:datetime(yyyy,12,31);
% % [si,ei]     = dateInds(datetime(yyyy,6,1),datetime(yyyy,9,1),tdaily);
% % vstart        = [vstart si];          
% % vcount        = [vcount (ei-si)];
% % 
% % 
% % % this was in compare_racmo_versions
% % 
% % % tried this to use ncrowcol before realizing it's not needed with my point
% % % interpolation approach
% % % get the lat/lon indices to read in the data around the catchment
% % % xpolymin    =   polymin.Vertices(:,1);
% % % xpolymax    =   polymax.Vertices(:,1);
% % % [slat,clat] =   ncrowcol(LAT,rot90(X,3),rot90(Y,3),xpolyb,ypolyb);
% % % [slon,clon] =   ncrowcol(LON,rot90(X,3),rot90(Y,3),xpolyb,ypolyb);
% % % [sswd,cswd] =   ncrowcol(swsd,rot90(X,3),rot90(Y,3),xpolyb,ypolyb);
% % 
% % % ... maybe I can just use the .nc files with the same workflow ...
% % 
% % 
% % 
% % 
% % 
% % 
% % % here's a way it could be done, from test_RACMO_refreeze
% % 
% % vers    = 'subsurface';
% % 
% % % set paths, get files
% % p.data  = ['/Users/coop558/mydata/racmo2.3/' vers '/'];
% % p.mask  = ['/Users/coop558/mydata/racmo2.3/' vers '/'];
% % list    = dir(fullfile([p.data '*.nc']));
% % 
% % load('projsipsn.mat');
% % load('racmo_subsurf_time.mat');
% % 
% % fnames      =   {list.name};
% % info        =   ncparse([p.data list(1).name]);
% % lon         =   readRacmo2p3([p.data list(1).name],'lon');
% % lat         =   readRacmo2p3([p.data list(1).name],'lat');
% % 
% % % is runoff ~= meltin + snowmelt - refreeze?
% % % get coordinates in the southwest
% % latdif      =   abs(lat-65);
% % londif      =   abs(lon+48);
% % [r,c]       =   find(latdif==min(latdif) & londif==min(londif));
% % t1          =   find(T==datetime(2016,7,1));
% % t2          =   find(T==datetime(2016,9,1));
% % tdaily           =   T(t1:t2);
% % ameltin     =   squeeze(meltin(r,c,t1:t2));
% % asnowmelt   =   squeeze(snowmelt(r,c,t1:t2));
% % arunoff     =   squeeze(runoff(r,c,t1:t2));
% % afreeze     =   squeeze(freeze(r,c,t1:t2));
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % % BELOW here is from the build metfiles script, above is testing to figure
% % % out how to do it directly from teh .nc file using ncrowcol
% % 
% % 
% % %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % %% paths
% % %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % p.data      =   '/Users/coop558/mydata/mar3.11/matfiles/region/level1/';
% % p.save      =   ['/Users/coop558/mydata/mar3.11/matfiles/region/level2/'...
% %                     albedo '/'];
% % p.modis     =   ['/Users/coop558/mydata/DATA/greenland/geus/albedo/'    ...
% %                     'region/matfiles/' proj '/native/'];
% % list        =   dir(fullfile([p.data '*.mat']));
% % 
% % %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                            
% % %% load the modis ice mask, extract x,y locations
% % %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % load([p.modis 'mask/modis_ice_mask']);
% % 
% % mask        =   icemask.mask;
% % elev        =   icemask.elev;
% % x           =   icemask.x;
% % y           =   icemask.y;
% % [nr,nc]     =   size(mask);
% % maskrs      =   logical(reshape(mask,nr*nc,1));
% % 
% % % reshape the elevation and x/y coord's and apply the mask
% % elev        =   reshape(elev,nr*nc,1);
% % x           =   reshape(x,nr*nc,1);
% % y           =   reshape(y,nr*nc,1);
% % x           =   x(maskrs);
% % y           =   y(maskrs);
% % elev        =   elev(maskrs);
% % %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % %% use this to figure out the indices for rio behar (it's 1045)
% % %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % % load('catchment_2016.mat');
% % % figure; 
% % % scatter(x,y,40,elev,'filled'); colorbar; hold on;
% % % plot(catchment.bestguess_x,catchment.bestguess_y);
% % % idx = find(x==-1.725e5 & y==-2.5125e6);
% % % [r,c] = ind2sub(size(icemask.mask),idx);
% % % scatter(x(idx),y(idx),400,'filled')
% % % plot(catchment.bestguess_x,catchment.bestguess_y);
% % %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % %% loop through all years and build metfiles for all the points
% % %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % for n = 1:nyears
% %     
% %     yyyy              =   yri+n-1;
% %     fname           =   [p.data 'mar_' int2str(yyyy)];
% %     load(fname)
% %     
% % %--------------------------------------------------------------------------
% %     % first read in modis albedo if requested and convert to hourly
% % %--------------------------------------------------------------------------    
% % 
% %   % build a daily calendar 
% %     startdate           =   [yyyy,1,1,12,0,0];
% %     enddate             =   [yyyy,12,31,12,0,0];
% %     t_daily             =   (datenum(startdate):1:datenum(enddate))';
% %     
% %   % build a 1-hourly calendar on to which the daily data will be resampled
% %     startdate           =   [yyyy,1,1,0,0,0];
% %     enddate             =   [yyyy,12,31,23,0,0];
% %     t_1hrly             =   (datenum(startdate):1/24:datenum(enddate))';
% %     nhrs                =   length(t_1hrly);
% % 
% %   % read in the modis data and interpolate it to 1-hourly
% %     if strcmp(albedo,'modis')
% %         load([p.modis 'modis_albedo_' int2str(yyyy)]);
% %         [nr,nc,nd]  =   size(modis.A);
% %         alb         =   reshape(modis.A,nr*nc,nd);
% %         alb         =   (alb(maskrs,:))';
% %         if ~isleap(yyyy)
% %             alb     =   alb(1:365,:);
% %         end
% %         alb         =   interp1(t_daily,alb,t_1hrly,'linear');
% %     else
% %         alb         =   (mar.albedo(maskrs,:))';
% %     end
% %     
% % %--------------------------------------------------------------------------
% %     % deal with the remaining variables
% % %--------------------------------------------------------------------------   
% % 
% %   % pull out the relevant forcing data and apply the mask 
% %     [Y,M,D,H,~,~]   =   datevec(mar.dates);
% %     swd             =   (mar.swd(maskrs,:))';
% %     lwd             =   (mar.lwd(maskrs,:))';
% %     tair            =   (mar.tair(maskrs,:))';
% %     relh            =   (mar.rh(maskrs,:))';
% %     wspd            =   (mar.wspd(maskrs,:))';
% %     wdir            =   (mar.wdir(maskrs,:))';
% %     ppt             =   (mar.prec(maskrs,:))';
% %     press           =   (mar.psfc(maskrs,:))';
% %     [nhrs,npts]     =   size(swd);
% % 
% %   % make the met files
% %     for n = 1:npts
% % 
% %         id_n        =   n.*ones(nhrs,1);
% %         x_n         =   x(n).*ones(nhrs,1);
% %         y_n         =   y(n).*ones(nhrs,1);
% %         elev_n      =   roundn(elev(n),0).*ones(nhrs,1);
% %         tair_n      =   linearfill(roundn(tair(:,n),-2));
% %         relh_n      =   linearfill(roundn(relh(:,n),0));
% %         wspd_n      =   linearfill(roundn(wspd(:,n),-1));
% %         wdir_n      =   linearfill(roundn(wdir(:,n),0));
% %         ppt_n       =   linearfill(roundn(ppt(:,n),-3));
% %         swd_n       =   linearfill(roundn(swd(:,n),-1));
% %         lwd_n       =   linearfill(roundn(lwd(:,n),-1));
% %         alb_n       =   linearfill(roundn(alb(:,n),-3));
% %         press_n     =   linearfill(roundn(press(:,n),-3));
% %         
% %         met         =   [Y,M,D,H,id_n,x_n,y_n,elev_n,tair_n,relh_n,     ...
% %                         wspd_n,wdir_n,ppt_n,swd_n,lwd_n,alb_n,press_n];
% % 
% %         if save_data == 1
% %             psave   =   [p.save int2str(yyyy)];
% %             if ~exist(psave,'dir')
% %                 mkdir(psave)
% %             end
% %             fsave   =   [psave '/met_' int2str(id_n(1)) '.mat'];
% %             save(fsave,'met')
% %         end
% %         clear met
% %         
% %     end
% %     
% % end
% % 
% % % make figures to confirm things make sense
% % % if plot_figs == 1
% % % figure
% % % scatter(x,y,20,swd(4000,:),'filled')
% % % title('SW Down'); colorbar
% % % 
% % % figure
% % % scatter(x,y,20,lwd(4000,:),'filled')
% % % title('LW Down'); colorbar
% % % 
% % % figure
% % % scatter(x,y,20,elev,'filled')
% % % title('Elevation'); colorbar
% % % end
% % 
% % 
% % 
% % % % I shoudl be able to get RH from this:
% % % Tf  =   273.16;
% % % A   =   6.1115 * 100.0;
% % % B   =   22.452;
% % % C   =   272.55;
% % % % Compute the saturated water vapor pressure at the surface.
% % % es0 =   A.*exp((B.*(mar.TTH - Tf))./(C + (mar.TTH - Tf))); 
% % % es  =   mar.QQH.*461.*mar.TTH;
% % % RH  =   es./es0;
