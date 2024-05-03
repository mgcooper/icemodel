function out = standardnames(whichdata)

   % Not all values have standard names, so I constructed some.

   % Define the grid and time dimension names.
   dimnames = {
      '', ... grid cell
      'latitude', ...
      'longitude', ...
      'projection_x_coordinate', ...
      'projection_y_coordinate', ...
      'height_above_mean_sea_level', ... 'altitude'
      'depth', ...
      'time'
      };

   % Define the variable names
   switch whichdata

      case 'ice1'

         % Note: I would need to redefine the data to conform to these:
         %
         % surface_temperature [K]
         % land_ice_surface_melt_flux [kg m2 s-1]
         % surface_snow_and_ice_refreezing_flux [kg m2 s-1]
         % land_ice_runoff_flux [kg m2 s-1]
         % surface_runoff_flux [kg m2 s-1]
         % subsurface_runoff_flux [kg m2 s-1]
         %
         % There are also "amounts" which are like my m w.e. but scaled by
         % density of water (mass per unit area):
         %
         % (there is no "land_ice_runoff_amount"
         % runoff_amount [kg m-2]
         % surface_runoff_amount [kg m-2]
         % subsurface_runoff_amount [kg m-2]
         %
         % There is also a standard name for the quantity:
         % surface_snow_and_ice_melt_flux
         % surface_snow_and_ice_refreezing_flux
         %
         % !! but no "land_ice_refreezing_flux"
         %
         % Note: this would define the ice surface temperature but not the snow
         % surface temperature if a snowpack were present (it is intended to
         % define the temperature that forces an ice sheet model):
         % temperature_at_top_of_ice_sheet_model [K]

         % Only surface_temperature is a recognized standard_name. The others
         % are modeled on similar official ones.
         out = {
            'surface_temperature', ...                               Tsfc
            '', ... 'cumulative_melt_amount_in_lwe', ...             melt
            '', ... 'cumulative_refreeze_amount_in_lwe', ...         freeze
            '', ... 'cumulative_sublimation_amount_in_lwe', ...      subl
            '', ... 'cumulative_condensation_amount_in_lwe', ...     cond
            '', ... 'cumulative_runoff_amount_in_lwe', ...           runoff
            '', ... 'cumulative_melt_amount_in_lwe', ...             depth_melt
            '', ... 'cumulative_refreeze_amount_in_lwe', ...         depth_freeze
            '', ... 'surface_runoff_amount_in_lwe', ...              surf_runoff
            '', ... 'subsurface_runoff_amount_in_lwe', ...           column_runoff
            };

      case 'ice2'

         % lwe_thickness_of_frozen_water_content_of_soil_layer [1]
         % volume_fraction_of_frozen_water_in_soil [1]
         % frozen_water_content_of_soil_layer [kg m-2]
         % land_ice_temperature [K]

         out = {
            'land_ice_temperature', ...
            '', ... 'change_in_volume_fraction_of_frozen_water_in_layer',...
            '', ... 'volume_fraction_of_frozen_water_in_layer', ...
            '', ... 'volume_fraction_of_unfrozen_water_in_layer', ...
            };

      case 'met'

      case {'dimensions', 'dims'}

         out = dimnames;

      otherwise
         error('unrecognized icemodel data file name')
   end

   % Return as a column
   out = out(:);
end
