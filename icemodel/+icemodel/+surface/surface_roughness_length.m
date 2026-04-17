function z0m = surface_roughness_length(snow_depth, ro_sfc, z0_ice, ...
      z0_snow_low_density, z0_snow_high_density)
   %SURFACE_ROUGHNESS_LENGTH Select the momentum roughness length z0m.
   %
   %  z0m = icemodel.surface.surface_roughness_length(snow_depth, ro_sfc, ...)
   %
   % Bare ice is selected by snow depth alone. The surface density only
   % splits the snow/firn branch into lower-density versus higher-density
   % roughness values.
   %
   %#codegen

   persistent ro_threshold
   if isempty(ro_threshold)
      ro_threshold = icemodel.parameterLookup( ...
         'thf_bulk_snow_density_threshold');
   end

   if snow_depth <= 0
      z0m = z0_ice;
      return
   end

   if ro_sfc <= ro_threshold
      z0m = z0_snow_low_density;
   else
      z0m = z0_snow_high_density;
   end
end
