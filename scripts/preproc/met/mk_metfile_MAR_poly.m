clean

% this script uses the new 'data' timetables to build catchment-scale
% metfiles. the only difference is that the metfile renames the variables
% for ease of use in the model, and there are options here to interpolate
% to higher resolution and save separate files with mar vs modis albedo 

% mk_metfile_MAR.m script calls makeMarMetfile.m function to build the
% point-scale metfiles. 


savedata =  true;
sitename =  'behar';
startyr  =  2015;
endyear  =  2016;
nyears   =  endyear-startyr+1;
interpT  =  true;
dtnew    =  '15m';

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% set paths
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pathdata =  ['/Users/coop558/mydata/mar3.11/matfiles/' sitename '/data/'];
pathsave =  setpath('GREENLAND/icemodel/input/met/');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% read in the catchment to get the x,y coordinate 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for n = 1:nyears
    
    thisyr  = num2str(startyr+n-1);
    fdata   = [pathdata 'MAR_' sitename '_' thisyr '.mat'];
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     met     = marData2Met(fdata);
    metcopy = met;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if savedata == true
        
        fsave = ['met_' sitename '_MAR_' thisyr];

        save([pathsave fsave '_1hr.mat'],'met');

        % interpolate, re-check out-of-bound values, and save
        if interpT == true
            met      = interpMet(met,dtnew);
            met.date = datenum(met.Time);
            met      = metchecks(met,false); % false = don't plot

            save([pathsave fsave '_' dtnew '.mat'],'met');
        end
        
% % JULY 2022, need to just add a column for MODIS or rely on the userData
% swap method, so i commented this out. 

%         % repeat, swapping out MAR albedo with MODIS albedo
%         met         = metcopy; clear metCopy
%         met.albedo  = met.MODIS;
%         
%         fsave = ['met_' sitename '_MAR_MODISalbedo_' thisyr];
% 
%         save([pathsave fsave '_1hr.mat'],'met');
%         save([p.copy fsave '_1hr.mat'],'met');
% 
%         % interpolate, re-check out-of-bound values, and save
%         if interpT == true
%             met      = interpMet(met,dtnew);
%             met.date = datenum(met.Time);
%             met      = metchecks(met,false); % false = don't plot
% 
%             save([pathsave fsave '_' dtnew '.mat'],'met');
% 
%         end        
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
% p.modis     =   ['/Users/coop558/mydata/DATA/greenland/geus/albedo/'    ...
%                     'region/matfiles/' proj '/native/'];
% list        =   dir(fullfile([pathdata '*.mat']));
% 
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                            
% %% load the modis ice mask, extract x,y locations
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% load([p.modis 'mask/modis_ice_mask']);
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
%         load([p.modis 'modis_albedo_' int2str(yyyy)]);
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
