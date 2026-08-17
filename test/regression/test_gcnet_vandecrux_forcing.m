function tests = test_gcnet_vandecrux_forcing
   %TEST_GCNET_VANDECRUX_FORCING Verify GC-Net/Vandecrux surface builders.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   % Install the canonical test path and allocate an isolated source cache.
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;
   testCase.TestData.tmp = tempname;
   mkdir(testCase.TestData.tmp)
end

function teardown(testCase)
   % Remove generated NetCDF fixtures.
   if isfolder(testCase.TestData.tmp)
      rmdir(testCase.TestData.tmp, 's')
   end
   clear testCase.TestData.cleanup
end

function test_buildGcnetVandecruxData_maps_units_and_metadata(testCase)
   % The Data builder maps native Vandecrux names, converts units, and records
   % the source-filled longwave/provenance policy on the userdata timetable.
   cache = makeGcnetSurfaceCache(testCase.TestData.tmp);

   [Data, metadata] = icemodel.forcing.buildGcnetVandecruxData( ...
      "dye2_long", source_dir=cache);

   expected_time = datetime(2000, 6, 1, 12, 0, 0, 'TimeZone', 'UTC') ...
      + hours(0:2)';
   testCase.verifyEqual(Data.Time, expected_time);
   testCase.verifyEqual(Data.tair, [260; 261; 262], 'AbsTol', 1e-12);
   testCase.verifyEqual(Data.psfc, [80000; 80100; 80200], 'AbsTol', 1e-12);
   testCase.verifyEqual(Data.snowf, [0; 0.001; 0.002] / 3600, ...
      'AbsTol', 1e-15);
   testCase.verifyEqual(Data.smb, [0.01; 0.02; 0.03], 'AbsTol', 1e-12);
   testCase.verifyEqual(Data.subl, [0.01; 0.02; 0.03], 'AbsTol', 1e-12);
   testCase.verifyEqual(Data.albedo, [0.7; 0.6; 0.5], 'AbsTol', 1e-12);
   testCase.verifyTrue(all(isnan(Data.rainf)));

   names = string(Data.Properties.VariableNames);
   units = string(Data.Properties.VariableUnits);
   testCase.verifyEqual(units(names == "snowf"), "m s-1");
   testCase.verifyEqual(units(names == "smb"), "mWE/h");
   testCase.verifyTrue(ismember("SRin_origin", metadata.source_variables));
   testCase.verifyEqual(size(metadata.source_variables, 2), 1);
   testCase.verifyTrue(all(isfinite([metadata.site_location.x_epsg3413, ...
      metadata.site_location.y_epsg3413])));
   testCase.verifyTrue(contains(metadata.lwd_policy, "source-filled"));
   testCase.verifyTrue(contains(metadata.lwd_policy, "HIRHAM5"));
   testCase.verifyTrue(contains(metadata.lwd_policy, "RACMO2.3p2"));
   testCase.verifyTrue(contains(metadata.precip_policy, ...
      "ppt = NaN placeholder"));
   testCase.verifyEqual(string(metadata.gcnet_subl_native_sign_convention), ...
      "negative_loss_positive_deposition");
   testCase.verifyEqual(string(metadata.gcnet_subl_sign_convention), ...
      "positive_loss_negative_deposition");
   testCase.verifyEqual(string(metadata.channel_map.lwd), "LRin");
   testCase.verifyEqual(string(Data.Properties.UserData.station), "DYE_2");
   testCase.verifyEqual(Data.Properties.CustomProperties.Elev, 2165);
   testCase.verifyEqual( ...
      metadata.albedo_transient_qc.episode_day_count, 0);
   testCase.verifyEmpty(metadata.albedo_transient_qc.episode_dates);
end

function test_buildGcnetVandecruxData_uses_normalized_cache_lookup(testCase)
   % Surface files accepted by fetch-style normalized matching must also build.
   cache = makeGcnetSurfaceCache(testCase.TestData.tmp);
   package_dir = fullfile(cache, 'surface-package');
   mkdir(package_dir);
   movefile(fullfile(cache, 'DYE_2_surface.nc'), ...
      fullfile(package_dir, 'dye-2_surface.nc'));

   Data = icemodel.forcing.buildGcnetVandecruxData( ...
      "dye2_long", source_dir=cache);

   testCase.verifyEqual(string(Data.Properties.UserData.station), "DYE_2");
   testCase.verifyEqual(Data.tair, [260; 261; 262], 'AbsTol', 1e-12);
end

function test_buildGcnetVandecruxData_rejects_implausible_accumulation_albedo(testCase)
   % Isolated ratios below the accumulation-site floor must not survive the
   % source adapter merely because downwelling shortwave exceeds 10 W m-2.
   cache = makeGcnetSurfaceCache(testCase.TestData.tmp);
   filename = fullfile(cache, 'DYE_2_surface.nc');
   ncwrite(filename, 'SRout', [20; 120; 150]);

   [Data, metadata] = icemodel.forcing.buildGcnetVandecruxData( ...
      "dye2_long", source_dir=cache);

   testCase.verifyTrue(isnan(Data.albedo(1)));
   testCase.verifyEqual(Data.albedo(2:3), [0.6; 0.5], 'AbsTol', 1e-12);
   testCase.verifyTrue(contains(metadata.albedo_policy, ...
      "albedo >= 0.3"));
   testCase.verifyEqual(metadata.albedo_qc_counts.below_minimum, 1);
   testCase.verifyEqual(metadata.albedo_qc_counts.total, 1);
end

function test_buildGcnetVandecruxData_uses_summit_dry_snow_floor(testCase)
   % Summit uses a regime-level dry-snow floor without censoring raw radiation.
   cache = makeGcnetSurfaceCache(testCase.TestData.tmp);

   [Data, metadata] = icemodel.forcing.buildGcnetVandecruxData( ...
      "summit", source_dir=cache);

   testCase.verifyEqual(Data.swu, [80; 130; 160], 'AbsTol', 1e-12);
   testCase.verifyEqual(Data.albedo(1), 80 / 110, 'AbsTol', 1e-12);
   testCase.verifyTrue(all(isnan(Data.albedo(2:3))));
   testCase.verifyTrue(contains(metadata.albedo_policy, ...
      "albedo >= 0.7"));
   testCase.verifyEqual(metadata.albedo_qc_counts.below_minimum, 2);
end

function test_buildGcnetVandecruxData_masks_transient_derived_episode(testCase)
   % A recovered daily reflected-shortwave collapse masks only derived channels;
   % the source swd/swu payload remains available for provenance and diagnosis.
   cache = makeGcnetTransientCache(testCase.TestData.tmp);

   [Data, metadata] = icemodel.forcing.buildGcnetVandecruxData( ...
      "summit", source_dir=cache);

   day = dateshift(Data.Time, 'start', 'day');
   expected_days = datetime(2000, 6, 1, 'TimeZone', 'UTC') ...
      + caldays([43; 45; 47]);
   flagged = ismember(day, expected_days);
   control = day == datetime(2000, 7, 15, 'TimeZone', 'UTC');
   testCase.verifyTrue(all(isfinite(Data.swd(flagged))));
   testCase.verifyTrue(all(isfinite(Data.swu(flagged))));
   testCase.verifyTrue(all(isnan(Data.albedo(flagged))));
   testCase.verifyTrue(all(isnan(Data.swn(flagged))));
   testCase.verifyTrue(all(isnan(Data.netr(flagged))));
   testCase.verifyTrue(any(isfinite(Data.albedo(control))));
   testCase.verifyEqual( ...
      metadata.albedo_transient_qc.seed_day_count, 1);
   testCase.verifyEqual( ...
      metadata.albedo_transient_qc.episode_day_count, 3);
   testCase.verifyEqual( ...
      string(metadata.albedo_transient_qc.episode_dates), ...
      string(expected_days, 'yyyy-MM-dd'));
   testCase.verifyTrue(contains(metadata.shortwave_balance_policy, ...
      "source swd/swu preserved"));
end

function test_gcnetTime_converts_units_and_rejects_unknown(testCase)
   % GC-Net builders and inventory use one time-unit parser.
   t = icemodel.forcing.helpers.gcnetTime( ...
      [0; 1], "hours since 2000-1-1 0:0:0");

   expected = datetime(2000, 1, 1, 0, 0, 0, 'TimeZone', 'UTC') ...
      + hours([0; 1]);
   testCase.verifyEqual(t, expected);
   testCase.verifyError( ...
      @() icemodel.forcing.helpers.gcnetTime(1, "seconds since 2000-1-1 0:0:0"), ...
      'icemodel:forcing:helpers:gcnetTime:unsupportedTimeUnits');
end

function test_buildGcnetMet_keeps_missing_precip_nan(testCase)
   % Native met has the required contract and canonical metadata, while ppt and
   % rainf stay missing because the source has snowfall but no rain channel.
   cache = makeGcnetSurfaceCache(testCase.TestData.tmp);

   [met, metadata, Data] = icemodel.forcing.buildGcnetVandecruxMet( ...
      "sum", source_dir=cache);

   required = icemodel.forcing.helpers.metvariables();
   names = string(met.Properties.VariableNames);
   testCase.verifyTrue(all(ismember(required, names)));
   testCase.verifyTrue(all(isnan(met.ppt)));
   testCase.verifyTrue(all(isnan(met.rainf)));
   testCase.verifyEqual(met.snowf, Data.snowf, 'AbsTol', 1e-15);
   testCase.verifyEqual(string(met.Properties.VariableUnits(names == "ppt")), ...
      "m s-1");
   testCase.verifyTrue(isprop(met.Properties.CustomProperties, ...
      "StandardNames"));
   testCase.verifyEqual( ...
      met.Properties.CustomProperties.StandardNames(names == "lwd"), ...
      "surface_downwelling_longwave_flux_in_air");
   testCase.verifyWarningFree(@() icemodel.forcing.helpers.validatemet(met));
   testCase.verifyTrue(contains(metadata.precip_policy, ...
      "rainf = NaN placeholder"));
   testCase.verifyTrue(contains( ...
      string(met.Properties.UserData.precip_policy), ...
      "ppt = NaN placeholder"));
   testCase.verifyTrue(metadata.fillwithmissing);
   testCase.verifyEqual(string(met.Properties.UserData.station), "Summit");
   testCase.verifyTrue(isfield(met.Properties.UserData, 'met_variables'));
end

function test_buildGcnetMet_supports_strict_mode(testCase)
   % fillwithmissing=false is exposed and leaves validation strict.
   cache = makeGcnetSurfaceCache(testCase.TestData.tmp);

   [met, metadata] = icemodel.forcing.buildGcnetVandecruxMet( ...
      "sum", source_dir=cache, fillwithmissing=false);

   testCase.verifyFalse(metadata.fillwithmissing);
   testCase.verifyWarningFree(@() icemodel.forcing.helpers.validatemet(met));
   testCase.verifyFalse(met.Properties.UserData.fillwithmissing);
end

function test_buildGcnetFirnTemp_maps_depth_and_provenance(testCase)
   % Firn-temperature observations are read as time x level canonical matrices
   % while retaining the original Vandecrux variable names and DOI provenance.
   cache = makeGcnetFirnTemperatureCache(testCase.TestData.tmp);

   [obs, metadata] = ...
      icemodel.forcing.buildGcnetVandecruxFirnTemperature( ...
      "dye2_long", source_dir=cache, ...
      startdate="2000-01-01 01:00:00", ...
      enddate="2000-01-01 02:00:00");

   expected_time = datetime(2000, 1, 1, 1, 0, 0, 'TimeZone', 'UTC') ...
      + hours(0:1)';
   testCase.verifyEqual(obs.time, expected_time);
   testCase.verifyEqual(obs.level, [1; 2]);
   testCase.verifySize(obs.subsurface_temperature, [2 2]);
   testCase.verifyEqual(obs.subsurface_temperature, ...
      [-9 -19; -8 -18] + 273.15, 'AbsTol', 1e-12);
   testCase.verifyEqual(obs.depth, [1.1 2.1; 1.2 2.2], ...
      'AbsTol', 1e-12);
   testCase.verifyEqual(string(obs.variables.subsurface_temperature), ...
      "T_firn");
   testCase.verifyEqual(string(obs.units.subsurface_temperature), "K");

   testCase.verifyEqual(string(metadata.station), "DYE_2");
   testCase.verifyTrue(any(string(metadata.aliases) == "dye2_long"));
   testCase.verifyEqual(string(metadata.provenance.doi), ...
      "10.18739/A2833N00P");
   testCase.verifyEqual(string(metadata.original_variables), ...
      ["T_firn"; "Depth"]);
   testCase.verifyTrue(contains(metadata.temperature_policy, "converted"));
   testCase.verifyEqual( ...
      string(metadata.source_variable_attributes.T_firn.units), "degC");
   testCase.verifyEqual(size(metadata.source_variables, 2), 1);
   testCase.verifyTrue(all(isfinite([metadata.site_location.x_epsg3413, ...
      metadata.site_location.y_epsg3413])));
end

function test_gcnet_forcing_name_sets_vandecrux_heights(testCase)
   % The staged gcnet met source uses the Vandecrux 2 m temperature/RH and
   % 10 m wind convention, not PROMICE single-boom heights.
   testCase.verifyTrue(ismember("gcnet", icemodel.namelists.forcings()));
   testCase.verifyWarningFree( ...
      @() icemodel.validators.mustBeForcingName("gcnet"));

   opts = icemodel.setopts('icemodel', 'dye2', 2000, 'gcnet');
   testCase.verifyEqual(opts.z_tair, 2.0);
   testCase.verifyEqual(opts.z_wind, 10.0);
   testCase.verifyEqual(opts.z_relh, opts.z_tair);
end

%% Local fixture helpers
function cache = makeGcnetSurfaceCache(root)
   %MAKEGCNETSURFACECACHE Create tiny Vandecrux surface NetCDF fixtures.
   cache = fullfile(root, 'gcnet-surface-cache');
   mkdir(cache);
   writeSurfaceFile(fullfile(cache, 'DYE_2_surface.nc'), 0);
   writeSurfaceFile(fullfile(cache, 'Summit_surface.nc'), 10);
end

function cache = makeGcnetTransientCache(root)
   %MAKEGCNETTRANSIENTCACHE Create a 90-day seeded-collapse source fixture.
   cache = fullfile(root, 'gcnet-transient-cache');
   mkdir(cache);
   writeTransientSurfaceFile(fullfile(cache, 'Summit_surface.nc'));
end

function cache = makeGcnetFirnTemperatureCache(root)
   %MAKEGCNETFIRNTEMPERATURECACHE Create a tiny firn-temperature fixture.
   cache = fullfile(root, 'gcnet-firn-temperature-cache');
   mkdir(cache);
   writeFirnTemperatureFile(fullfile(cache, 'DYE_2_T_firn_obs.nc'));
   writeFirnTemperatureXml(cache);
end

function writeSurfaceFile(filename, offset)
   %WRITESURFACEFILE Write the source variables used by the GC-Net builders.
   time = 36676 + 0.5 + (0:2) / 24;
   writeVar(filename, "time", time, "days since 1900-1-1 0:0:0", ...
      "Time (UTC)");

   writeVar(filename, "SRin", [100 200 300] + offset, "Wm-2", ...
      "Incoming shortwave radiation");
   writeVar(filename, "SRin_origin", [0 101 0], "-", "Origin of SRin");
   writeVar(filename, "SRout", [70 120 150] + offset, "Wm-2", ...
      "Outgoing shortwave radiation");
   writeVar(filename, "SRout_origin", [0 101 0], "-", "Origin of SRout");
   writeVar(filename, "LRin", [250 251 252] + offset, "Wm-2", ...
      "Downward longwave radiation flux");
   writeVar(filename, "LRout", [240 241 242] + offset, "Wm-2", ...
      "Upward longwave radiation flux");
   writeVar(filename, "SHF", [1 2 3], "Wm-2", "Sensible heat flux");
   writeVar(filename, "LHF", [4 5 6], "Wm-2", "Latent heat flux");
   writeVar(filename, "H_surf_obs", [10 10.1 10.2], "m", ...
      "Observed surface height");
   writeVar(filename, "SMB", [0.01 0.02 0.03], "m_weq", ...
      "Surface mass balance");
   writeVar(filename, "melt", [0.1 0.2 0.3], "m_weq", "Melt");
   writeVar(filename, "sublimation", [-0.01 -0.02 -0.03], "m_weq", ...
      "Net sublimation");
   writeVar(filename, "snowfall", [0 0.001 0.002], "m_weq", "Snowfall");
   writeVar(filename, "snowfall_origin", [0 101 101], "-", ...
      "Origin of snowfall estimate");
   writeVar(filename, "Tsurf", [259 260 261], "K", "Surface temperature");
   writeVar(filename, "Ta_2m", [260 261 262], "K", "2m air temperature");
   writeVar(filename, "Ta_2m_origin", [0 0 101], "-", "Origin of Ta");
   writeVar(filename, "RH_2m", [70 71 72], "%", ...
      "2m relative humidity with regards to water");
   writeVar(filename, "RH_2m_origin", [0 0 101], "-", "Origin of RH");
   writeVar(filename, "WS_10m", [5 6 7], "ms-1", "10 m wind speed");
   writeVar(filename, "WS_10m_origin", [0 0 101], "-", "Origin of WS");
   writeVar(filename, "Pres", [800 801 802], "hPa", ...
      "air pressure at the surface");
end

function writeTransientSurfaceFile(filename)
   %WRITETRANSIENTSURFACEFILE Write daily context around one collapse episode.
   n = 90 * 24;
   time = 36676 + (0:n - 1) / 24;
   writeVar(filename, "time", time, "days since 1900-1-1 0:0:0", ...
      "Time (UTC)");

   % Constant irradiance isolates the cross-day albedo rule from solar geometry.
   swd = 100 * ones(1, n);
   alpha = 0.85 * ones(1, n);
   day_index = floor((0:n - 1) / 24) + 1;
   alpha(day_index == 44) = 0.72;
   alpha(day_index == 46) = 0.40;
   alpha(day_index == 48) = 0.72;
   writeVar(filename, "SRin", swd, "Wm-2", ...
      "Incoming shortwave radiation");
   writeVar(filename, "SRout", alpha .* swd, "Wm-2", ...
      "Outgoing shortwave radiation");
   writeVar(filename, "LRin", 250 * ones(1, n), "Wm-2", ...
      "Downward longwave radiation flux");
   writeVar(filename, "LRout", 240 * ones(1, n), "Wm-2", ...
      "Upward longwave radiation flux");
end

function writeVar(filename, name, values, units, long_name)
   %WRITEVAR Write one one-dimensional NetCDF fixture variable.
   nccreate(filename, name, 'Dimensions', {'time', numel(values)});
   ncwrite(filename, name, values);
   ncwriteatt(filename, name, 'units', units);
   ncwriteatt(filename, name, 'long_name', long_name);
end

function writeFirnTemperatureFile(filename)
   %WRITEFIRNTEMPERATUREFILE Write a level x time firn-temperature fixture.
   time = [876576 876577 876578];
   nccreate(filename, "time", 'Dimensions', {'time', numel(time)});
   ncwrite(filename, "time", time);
   ncwriteatt(filename, "time", 'units', "hours since 1900-1-1 0:0:0");

   nccreate(filename, "level", 'Dimensions', {'level', 2});
   ncwrite(filename, "level", [1 2]);

   temp = [-10 -9 -8; -20 -19 -18];
   nccreate(filename, "T_firn", 'Dimensions', {'level', 2, 'time', 3});
   ncwrite(filename, "T_firn", temp);
   ncwriteatt(filename, "T_firn", 'units', "degC");
   ncwriteatt(filename, "T_firn", 'long_name', "firn temperature");

   depth = [1 1.1 1.2; 2 2.1 2.2];
   nccreate(filename, "Depth", 'Dimensions', {'level', 2, 'time', 3});
   ncwrite(filename, "Depth", depth);
   ncwriteatt(filename, "Depth", 'units', "m");
   ncwriteatt(filename, "Depth", 'long_name', ...
      "measurement depth below the surface");
end

function writeFirnTemperatureXml(cache)
   %WRITEFIRNTEMPERATUREXML Write the minimal Arctic Data DOI metadata fixture.
   text = [ ...
      '<eml packageId="doi:10.18739/A2833N00P">' ...
      '<dataset>' ...
      '<title>Firn temperatures and measurement depths at nine sites</title>' ...
      '<abstract><para>Vandecrux et al. 2020.</para></abstract>' ...
      '</dataset>' ...
      '</eml>'];
   fid = fopen(fullfile(cache, ...
      'Firn_temperatures_and_measurement_depths_at_nine.xml'), 'w');
   fwrite(fid, text);
   fclose(fid);
end
