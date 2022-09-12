clean

%~~~~~~~~~~~~~~~~~~~~~~~~~~~  set options
save_data   =   true;
y1          =   2009;
y2          =   2018;
basin       =   'ak4';
models      =   {'skinmodel','skinmodel','icemodel','icemodel'};
albedos     =   {'mar','modis','mar','modis'};
climate     =   'reference';
res         =   250;                    % resampling resolution
buffer      =   7500;                   % this requires guess and check. 
testbuffer  =   false;                  % set false once buffer is correct 
method      =   'natural';
nyears      =   y2-y1+1;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~  set paths
p.shape     = ['/Volumes/Samsung_T5/DATA/greenland/runoff/' basin '/'];
p.data      = ['/Volumes/Samsung_T5/matlab/GREENLAND/runoff/icemodel/' ...
                'data/output/region/'];
p.save      = ['/Volumes/Samsung_T5/matlab/GREENLAND/runoff/icemodel/' ...
                'data/output/' basin '/'];

%~~~~~~~~~~~~~~~~~~~~~~~~~~~ get the catchment boundary
load('projsipsn.mat')
f       =   'd8_basin_ice_vec_kbAK4_fw1p0';
f       =   [p.shape f '/d8_basin_vec_kbAK4_fw1p0.shp'];
sf      =   get_catchment(f,projsipsn,'no plot');

%~~~~~~~~~~~~~~~~~~~~~~~~~~ get the rcm grid points from the ice mask 
load('modis_ice_mask'); 
x       =   icemask.x(logical(icemask.mask));
y       =   icemask.y(logical(icemask.mask));

%~~~~~~~~~~~~~~~~~~~~~~~~~~ get interpolation points within the catchment
poly    =   sf.poly;
pts     =   interpolation_points(poly,x,y,buffer,res,'no plot');
mask    =   pts.maskbuffer.inpolyb;
x       =   pts.xinbuffer;              % rcm grid points
y       =   pts.yinbuffer;
xq      =   pts.xq;                     % interpolation grid points
yq      =   pts.yq;

%~~~~~~~~~~~~~~~~~~~~~~~~~~ interpolate the rcm variables to the fine grid

if testbuffer == false
    
for i = 1:length(models)
    
    model   =   models{i};
    albedo  =   albedos{i};
    
    for j = 1:nyears

        % load the icemodel data
        yy      =   int2str(y1+j-1);
        f       =   [p.data model '_' climate '_' albedo '_albedo_' yy]; 
        load(f)
        
        % clip out the rcm grid points using the mask
        r       =   data.R(:,mask)';
        T       =   data.Time; clear data

        % interpolate the rcm data to the fine grid
        R        =   interpolate_rcm(x,y,r,xq,yq,method);
        R.Time   =   T;
        R.units  =   'mWE';

        % save the data
        if save_data == true
            f   = [p.save model '_' climate '_' albedo '_albedo_' yy];
            save(f,'R');
        end
        
        clear R r
    end    
end

end




