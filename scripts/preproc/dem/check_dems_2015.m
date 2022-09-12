% catch all script to process snowmodel input for Rio Behar

% this script looks at the 'SPIRIT' dems in Kang's folder. The resolution
% is 40 m. The GIMP dem I used to run snowmodel back in 2015 has 30 m
% resolution. I might as well just use that dem for all years for the
% purposes of SnowModel, to keep things equal. 

clean

%%
path.dem            =   ['/Users/mattcooper/Google UCLA/temp/Greenland/' ...
                            'rio_behar/dem/2015/'];
                
%%

% 2015 study as 20-23 June

%%
list                =   dir(fullfile([path.dem '*.tif']));
nodata              =   65535;

for n = 1:length(list)
    
    % get the info and filename
    info{n}         =   geotiffinfo([path.dem list(n).name]);
    f_n             =   [list(n).name(9:end-4) '_' list(n).name(1:7)];
    
    % build unique names for the sturcutre
    fdem_n          =   ['dem' int2str(n)];
    fR_n            =   ['R' int2str(n)];
    fname_n         =   ['fname' int2str(n)];
    
    % get the data
    [dem_n,R_n]     =   geotiffread([path.dem list(n).name]);
    dem_n           =   double(dem_n);
    bi_n            =   dem_n == nodata;
    
    % put the data in the structure
    dem_n(bi_n)     =   nan;
    dem.(fdem_n)    =   dem_n;
    dem.(fR_n)      =   R_n;
    dem.(fname_n)   =   f_n;
end

%%

for n = 1:length(list)
    fdem_n          =   ['dem' int2str(n)];
    fR_n            =   ['R' int2str(n)];
    f(n)            =   figure;
    m(n)            =   rastersurf(dem.(fdem_n),dem.(fR_n));
end

%% 

% dem 4 is SPIRIT_DEM_WQ7_ps_08JUN27. I want to use this one
% the resolution is 40 m, which is plenty





