clean

% There was a version of this script called _test that tested only writing the
% 2d vars (and maybe only the 0-d vars based on the note below but the note may
% be a mistake, it seems it tested only writing 2-d vars). In the test version,
% comprsz was set to 9, so I added that to the if test_compression == true
% section, but I think that was just a one-off setting, separate from the actual
% test of including 0-d, 1-d, or 2-d vars.

savedata = true;
sitename = 'sector';
siteopts = setBasinOpts('sitename', sitename, 'simmodel', 'icemodel');
modelver = 'v10b';
simmodel = {'icemodel', 'skinmodel'};
userdata = {'mar', 'modis'};

% I also tested if compression is better w/o the 1-d and 2-d vars - it isn't
test_compression = false;
if test_compression == true
   comprsz = 9;           % compression size
else
   comprsz = 1;           % compression size
end
% file size in mb/cell for min/med/max compression
% 1: 10.9 (10.8 w/o 1-d or 0-d vars)
% 5: 9.4
% 9: 9.2 (9.1 w/o 1-d or 0-d vars)

simyears = siteopts.simyears;
numyears = numel(simyears);

% it's best to set these here rather than programatically
nvars0d = 6;  % lat, lon, x, y, elev, time
nvars1d = 4;  % Tsfc, melt, freeze, runoff
nvars2d = 4;  % Tice, f_liq, f_ice, df_liq

%% Set path to data
pathdata = fullfile(getenv('ICEMODELOUTPUTPATH'), sitename);
pathsave = fullfile(getenv('ICEMODELOUTPUTPATH'), sitename, 'nc');
pathmet  = fullfile(getenv('ICEMODELINPUTPATH'), 'met', sitename);

if ~isfolder(pathsave)
   mkdir(pathsave);
end

% set the dimesions
ncells = 3; % 1487    % store the data as a list
nlyrs = 500;         % number of vertical layers
nhrs = 8760;        % number of hours per year

% set default format
netcdf.setDefaultFormat('NC_FORMAT_NETCDF4');

%% Define variables, dimensions, put the data, close the file

if test_compression == true

   varnames = {'f_ice','f_liq','df_liq','Tice'};
   longnames = {'fraction of frozen water in control volume',            ...
      'fraction of unfrozen water in control volume',          ...
      'change in fraction of unfrozen water in control volume',...
      'thermodynamic temperature of control volume'};
   units = {'1','1','1','K'};
else

   varnames = {'f_ice','f_liq','df_liq','Tice','Tsfc','depth_melt',     ...
      'depth_freeze','column_runoff','latitude','longitude',   ...
      'x_easting','y_northing','elevation','time'};
   longnames = {'fraction of frozen water in control volume',            ...
      'fraction of unfrozen water in control volume',          ...
      'change in fraction of unfrozen water in control volume',...
      'thermodynamic temperature of control volume',           ...
      'thermodynamic temperature of ice surface',              ...
      'cumulative melt',                                       ...
      'cumulative refreeze',                                   ...
      'cumulative runoff',                                     ...
      'latitude of grid cell center',                          ...
      'longitude of grid cell center',                         ...
      'x-coordinate of grid cell center in EPSG:3413 WGS 84 / NSIDC Sea Ice Polar Stereographic North', ...
      'y-coordinate of grid cell center in EPSG:3413 WGS 84 / NSIDC Sea Ice Polar Stereographic North', ...
      'elevation above mean sea level',                        ...
      'seconds since 1 January, 00:00'};
   units = {'1','1','1','K','K','m w.e.','m w.e.','m w.e.',         ...
      'degree_north','degree_east','m','m','m','s'};
end

% should be: volume fraction of frozen water, volume fraction of unfrozen water
% standardnames = {'volume fraction of frozen water','','','land_ice_temperature'};

nvars = numel(varnames);

% CREATE THE FILES

for m = 1:1 %nyrs

   thisyear = num2str(simyears(m));

   filename = fullfile(pathsave, ['icemodel_' thisyear '.nc4']);

   % load one ice1 file to get the calendar
   load(fullfile(pathdata, thisyear, 'ice1_1.mat'), 'ice1');

   % create the file and set 'NOFILL' to improve performance
   ncid = netcdf.create(filename,'NETCDF4');
   netcdf.setFill(ncid,'NC_NOFILL');

   % define the dimensions. the unlimited dimension doesn't seem to be
   % necessary, i got a few errors that appeared to be fixed by using it,
   % but then once the script was finished, it works both ways. the unidata
   % blog says compression may be degraded if unlimited is used, but the
   % file sizes in my tests were identical

   % dimidcells = netcdf.defDim(ncid,'cells', netcdf.getConstant('NC_UNLIMITED'));
   dimidcells = netcdf.defDim(ncid,'cells', ncells);
   dimidlayers = netcdf.defDim(ncid,'layers',nlyrs);
   dimidtime = netcdf.defDim(ncid,'time',nhrs);

   dimid0d = dimidcells;
   dimid1d = [dimidcells dimidtime];
   dimid2d = [dimidcells dimidlayers dimidtime];

   if test_compression == true

      dimids = {dimid2d, dimid2d, dimid2d, dimid2d};
   else
      dimids = {dimid2d, dimid2d, dimid2d, dimid2d, dimid1d, dimid1d, ...
         dimid1d, dimid1d, dimid0d, dimid0d, dimid0d, dimid0d, ...
         dimid0d, dimidtime};
   end

   % ALL VARS - set the varid's and attributes
   varid = nan(nvars,1);

   if test_compression == true
      writevars_test(pathdata, thisyear, ...
         varid, nvars, varnames, ncid, dimids, comprsz, ...
         longnames, units)
   else

      for p = 1:nvars
   
         thisvar = varnames{p};
         varid(p) = netcdf.defVar(ncid,thisvar,'NC_DOUBLE',dimids{p});
   
         netcdf.defVarDeflate(ncid,varid(p),true,true,comprsz);
         netcdf.putAtt(ncid,varid(p),'long_name',longnames{p});
         netcdf.putAtt(ncid,varid(p),'units',units{p});
      end
   
   
      for p = 1:nvars2d + nvars1d
   
         thisvar = varnames{p};
         thisid = varid(p);
   
         for n = 1:ncells
   
            load(fullfile(pathdata, thisyear, ['ice1_' num2str(n) '.mat']), 'ice1');
            load(fullfile(pathdata, thisyear, ['ice2_' num2str(n) '.mat']), 'ice2');
   
            % for 2d vars, need to transpose and set nlyrs start,count
            if p <= nvars2d
               data = transpose(rmleapinds(ice2.(thisvar),ice1.Time));
               netcdf.putVar(ncid,thisid,[n-1 0 0],[1 nlyrs nhrs],data);
            else
               data = rmleapinds(ice1.(thisvar),ice1.Time);
               netcdf.putVar(ncid,thisid,[n-1 0],[1 nhrs],data);
            end
         end % ncells
      end % nvars
   end
end % nyears

%% 0-d variables

if test_compression == true
   % do nothing, only write 2-d vars.
else

   % init the 0-d arrays
   vars0d.latitude = nan(ncells,1);
   vars0d.longitude = nan(ncells,1);
   vars0d.x_easting = nan(ncells,1);
   vars0d.y_northing = nan(ncells,1);
   vars0d.elevation = nan(ncells,1);
   vars0d.time = 0:3600:3600*nhrs-1;

   % read in the metfiles to get the lat/lon/x/y/elev
   for n = 1:ncells


      metfile = ['met_' num2str(n) '.mat'];
      load(fullfile(pathmet, metfile), 'met');

      vars0d.latitude(n) = met.Properties.CustomProperties.Lat;
      vars0d.longitude(n) = met.Properties.CustomProperties.Lon;
      vars0d.x_easting(n) = met.Properties.CustomProperties.X;
      vars0d.y_northing(n) = met.Properties.CustomProperties.Y;
      vars0d.elev(n) = met.Properties.CustomProperties.Elev;
   end

   % write the lat, lon, time vars
   for n = 1:nvars0d-1
      thisid = nvars1d+nvars2d+n;
      thisvar = varnames{thisid};
      netcdf.putVar(ncid,varid(thisid),vars0d.(thisvar));
   end
end

netcdf.close(ncid);

%% Check the result

info = ncinfo(filename);
test = ncread(filename,'f_liq');
test1 = squeeze(test(1,:,:));
test2 = squeeze(test(2,:,:));
test3 = squeeze(test(3,:,:));

figure; plot(test1(1,:)); hold on; plot(test1(2,:)); plot(test1(3,:));


% % % % % % % % % % % % % % % % % % % % % % % % % %
% % % this shows that runoff cannot be computed from f_liq after the fact,
% % whereas it can be computed from df_liq, so I need to save df_liq. I
% % should be able to compute ro_sno, cp_sno, etc from f_liq and f_ice, so I
% % still need to save them along with Tice
%
%
% vzero = zeros(size(ice2.f_liq(:,1)));
% df_test = [vzero diff(ice2.f_liq,1,2)];
%
% dz = 0.04;
% runoff_test = tocolumn(cumsum( sum( 4*dz.*df_test ) ));
% % runoff_test = tocolumn(cumsum( sum( dz.*ice2.df_liq ) ));
%
% figure;
% plot(ice1.Time,ice1.runoff); hold on;
% plot(ice1.Time,4.*runoff_test,':');

function writevars_test(pathdata, thisyear, varid, nvars, varnames, ...
      ncid, dimids, comprsz, longnames, units)

   for p = 1:nvars

      thisvar = varnames{p};
      varid(p) = netcdf.defVar(ncid,thisvar,'NC_DOUBLE',dimids{p});

      netcdf.defVarDeflate(ncid,varid(p),true,true,comprsz);
      netcdf.putAtt(ncid,varid(p),'long_name',longnames{p});
      netcdf.putAtt(ncid,varid(p),'units',units{p});

      for n = 1:ncells

         load(fullfile(pathdata, thisyear, ['ice1_' num2str(n) '.mat']), 'ice1');
         load(fullfile(pathdata, thisyear, ['ice2_' num2str(n) '.mat']), 'ice2');

         % for 2d vars, need to transpose and set nlyrs start,count
         data = transpose(rmleapinds(ice2.(thisvar),ice1.Time));
         netcdf.putVar(ncid,varid(p),[n-1 0 0],[1 nlyrs nhrs],data);
      end % ncells
   end % nvars
end
