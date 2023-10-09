clean

savedata    =  false;
sitename    =  'region';
yri         =  2008;
yrf         =  2018;
nyears      =  yrf-yri+1;
tinterp     =  true;
newdt       =  '15m';

%------------------------------------------------------------------------------
% set paths
%------------------------------------------------------------------------------

pathdata    =  setpath('mar3.11','data');
pathsave    =  setpath('mar3.11/matfiles/region/level2/grids/','data');
pathmodis   =  setpath('geus/albedo/gris/','data');
filelist    =  getlist(pathdata,'*.nc');

%------------------------------------------------------------------------------
% read in the catchment to get the x,y coordinate 
%------------------------------------------------------------------------------

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
% n=nyears;

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

%------------------------------------------------------------------------------
        
%    [met,data]  =  makeMarMetfileGrid(fmar,'Xmar',Xmar,'Ymar',Ymar,      ...
%                   'ROI',ROI,'Xmet',Xmet,'Ymet',Ymet,'fmodis',fmodis,    ...
%                   'Xmodis',Xmod,'Ymodis',Ymod);
               
% for building at the native mar grid cell resolution:
   [met,data]  =  makeMarMetfileGrid(fmar,'Xmar',Xmar,'Ymar',Ymar,      ...
                  'ROI',ROI,'fmodis',fmodis,'Xmodis',Xmod,'Ymodis',Ymod);

%------------------------------------------------------------------------------

   if savedata == true
      save(fsavemet,'met','-v7.3')
      save(fsavedata,'data','-v7.3')
   end

   clear met data
   
end
