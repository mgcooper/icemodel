function map = variables()
   %VARIABLES Canonical variable-metadata map for every icemodel channel.
   %
   %  map = icemodel.netcdf.defaults.variables()
   %
   % Returns a dictionary keyed by convenience name (the channel name the
   % forcing builders, met files, and netcdf writers use) whose values are
   % structs with the fields:
   %
   %    standard_name : CF standard_name (validated against the official CF
   %                    Standard Name Table) or '' when no CF name applies
   %    long_name     : human-readable description (CF long_name)
   %    unit          : canonical unit string
   %    is_cf         : true when standard_name is in the official CF table
   %
   % This is the ONE source of truth that replaces the per-dataset hand maps
   % in icemodel.netcdf.defaults.{standardnames,longnames,units} and the
   % duplicated unit list in icemodel.forcing.helpers.variableUnits, which is
   % now a thin wrapper over this map. Access a single channel with
   % icemodel.netcdf.defaults.variable(name).
   %
   % CF standard names here are validated programmatically against the
   % official table (icemodel.netcdf.defaults.cfStandardNames); see the
   % validatecf=true path of icemodel.netcdf.defaults.variable. Where a
   % channel has no official CF name (model diagnostics, instrument
   % channels), standard_name is '' and is_cf is false.
   %
   % See also: icemodel.netcdf.defaults.variable,
   %  icemodel.netcdf.defaults.cfStandardNames,
   %  icemodel.forcing.helpers.variableUnits

   % name, standard_name, long_name, unit
   % standard_name == '' marks a non-CF (hand-constructed or absent) name.
   rows = {
   % --- forcing / met channels --------------------------------------------
   'tair',   'air_temperature',                            'near-surface air temperature',                       'K'
   'tsfc',   'surface_temperature',                        'surface skin temperature',                           'K'
   'swd',    'surface_downwelling_shortwave_flux_in_air',  'downwelling shortwave radiative flux',               'W m-2'
   'swu',    'surface_upwelling_shortwave_flux_in_air',    'upwelling shortwave radiative flux',                 'W m-2'
   'lwd',    'surface_downwelling_longwave_flux_in_air',   'downwelling longwave radiative flux',                'W m-2'
   'lwu',    'surface_upwelling_longwave_flux_in_air',     'upwelling longwave radiative flux',                  'W m-2'
   'swn',    'surface_net_downward_shortwave_flux',        'net downward shortwave radiative flux',              'W m-2'
   'lwn',    'surface_net_downward_longwave_flux',         'net downward longwave radiative flux',               'W m-2'
   'netr',   'surface_net_downward_radiative_flux',        'net downward radiative flux',                        'W m-2'
   'shf',    'surface_downward_sensible_heat_flux',        'surface downward sensible heat flux',                'W m-2'
   'lhf',    'surface_downward_latent_heat_flux',          'surface downward latent heat flux',                  'W m-2'
   'thf',    'surface_downward_heat_flux_in_air',          'surface downward turbulent heat flux (sensible+latent)', 'W m-2'
   'albedo', 'surface_albedo',                             'surface broadband albedo',                           '-'
   'modis',  '',                                           'MODIS surface broadband albedo',                     '-'
   'cfrac',  'cloud_area_fraction',                        'cloud area fraction',                                '-'
   'rh',     'relative_humidity',                          'relative humidity',                                  '%'
   'rh_ice', '',                                           'relative humidity with respect to ice (subfreezing) or water', '%'
   'wspd',   'wind_speed',                                 'wind speed',                                         'm s-1'
   'wdir',   'wind_from_direction',                        'wind direction (from, clockwise from north)',        'degrees'
   'uwind',  'eastward_wind',                              'eastward wind component',                            'm s-1'
   'vwind',  'northward_wind',                             'northward wind component',                           'm s-1'
   'psfc',   'surface_air_pressure',                       'surface air pressure',                               'Pa'
   'shum',   'specific_humidity',                          'specific humidity',                                  'kg/kg'
   % --- precipitation (water-equivalent rate, m s-1) ----------------------
   'ppt',    'precipitation_flux',                         'total precipitation as water-equivalent rate',       'm s-1'
   'precip', 'precipitation_flux',                         'total precipitation as water-equivalent rate',       'm s-1'
   'rainf',  'rainfall_flux',                              'rainfall as water-equivalent rate',                  'm s-1'
   'snowf',  'snowfall_flux',                              'snowfall as water-equivalent rate',                  'm s-1'
   'rain',   'rainfall_flux',                              'rainfall as water-equivalent rate',                  'm s-1'
   'snow',   'snowfall_flux',                              'snowfall as water-equivalent rate',                  'm s-1'
   % --- mass-flux diagnostics (mWE/h) -------------------------------------
   'melt',     'surface_snow_and_ice_melt_flux',           'surface snow and ice melt flux',                     'mWE/h'
   'runoff',   'surface_runoff_flux',                      'surface runoff flux',                                'mWE/h'
   'evap',     'water_evaporation_flux',                   'water evaporation flux',                             'mWE/h'
   'smb',      'land_ice_surface_specific_mass_balance_rate', 'land-ice surface specific mass balance rate',     'mWE/h'
   'refreeze', 'surface_snow_and_ice_refreezing_flux',     'surface snow and ice refreezing flux',               'mWE/h'
   'subl',     'surface_snow_and_ice_sublimation_flux',    'surface snow and ice sublimation flux',              'mWE/h'
   'sndiv',    '',                                          'snow mass divergence (blowing-snow transport)',      'mWE/h'
   'meltin',   '',                                          'meltwater input to the snowpack',                    'mWE/h'
   % --- stores / heights / depths -----------------------------------------
   'swe',              '', 'snow water equivalent',                           'kg m-2'
   'snowd',            'surface_snow_thickness', 'snow depth',                'm'
   'snow_depth',       'surface_snow_thickness', 'snow depth',                'm'
   'ablation',         '', 'surface ablation (lowering) height',              'm'
   'surface_height',   '', 'net surface-height change relative to window start (accumulation sites; positive up)', 'm'
   'surface_height_flag', '', 'quality flag for the surface-height channel (0=direct observation, 1=gap-bridged/interpolated)', '1'
   'station_transition_flag', '', 'quality flag marking samples within a station-handover (AWS merge) window (0=no, 1=yes)', '1'
   'step_detected_flag', '', 'quality flag marking a detected single-timestep step-shift candidate (0=no, 1=yes)', '1'
   'step_correctable_flag', '', 'quality flag marking an UNAMBIGUOUS (correctable) step-shift (0=no, 1=yes)', '1'
   'step_magnitude',   '', 'signed magnitude of a detected surface-series step-shift at this sample (0 elsewhere)', 'm'
   'boom_height',      '', 'AWS boom height above the surface',               'm'
   'stake_height',     '', 'ablation-stake height above the surface',         'm'
   'transducer_depth', '', 'sonic-ranger / transducer depth to the surface',  'm'
   'elev',             'height_above_mean_sea_level', 'elevation above mean sea level', 'm'
   % --- corrected shortwave (PROMICE tilt/bias) ---------------------------
   'swd_cor', 'surface_downwelling_shortwave_flux_in_air', 'tilt/bias-corrected downwelling shortwave flux', 'W m-2'
   'swu_cor', 'surface_upwelling_shortwave_flux_in_air',   'tilt/bias-corrected upwelling shortwave flux',   'W m-2'
   % --- instrument geometry (PROMICE) -------------------------------------
   'tilt_x', '', 'platform tilt about the x axis',  'degrees'
   'tilt_y', '', 'platform tilt about the y axis',  'degrees'
   'rot',    '', 'platform azimuth (rotation)',     'degrees'
   % --- firn / ice-temperature string -------------------------------------
   'tice10m', 'land_ice_temperature', '10 m subsurface (firn) temperature', 'K'
   % --- ice1 model-output diagnostics (m w.e. cumulative) -----------------
   'melt_cum',      '', 'cumulative melt amount in meters liquid water equivalent',                'm w.e.'
   'freeze',        '', 'cumulative refreeze amount in meters liquid water equivalent',            'm w.e.'
   'cond',          '', 'cumulative condensation amount in meters liquid water equivalent',        'm w.e.'
   'runoff_cum',    '', 'cumulative runoff amount in meters liquid water equivalent',              'm w.e.'
   'depth_melt',    '', 'cumulative melt amount in meters liquid water equivalent',                'm w.e.'
   'depth_freeze',  '', 'cumulative refreeze amount in meters liquid water equivalent',            'm w.e.'
   'surf_runoff',   '', 'cumulative surface runoff amount in meters liquid water equivalent',      'm w.e.'
   'column_runoff', '', 'cumulative subsurface runoff amount in meters liquid water equivalent',   'm w.e.'
   % --- ice2 model-output diagnostics -------------------------------------
   'Tice',  'land_ice_temperature', 'thermodynamic temperature of control volume (layer)', 'degC'
   'f_ice', '', 'volume fraction of frozen water in control volume (layer)',   '1'
   'f_liq', '', 'volume fraction of unfrozen water in control volume (layer)', '1'
   'df_liq', '', 'change in volume fraction of unfrozen water in control volume (layer)', '1'
   };

   names = string(rows(:, 1));
   map = dictionary();
   for k = 1:size(rows, 1)
      s.standard_name = char(rows{k, 2});
      s.long_name = char(rows{k, 3});
      s.unit = char(rows{k, 4});
      s.is_cf = ~isempty(s.standard_name);
      map(names(k)) = s;
   end
end
