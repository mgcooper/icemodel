function ro_sfc = surface_bulk_density(f_ice_sfc, f_liq_sfc)
   %SURFACE_BULK_DENSITY Compute the bulk density of the top model layer.
   %
   %  ro_sfc = icemodel.surface.surface_bulk_density(f_ice_sfc, f_liq_sfc)
   %
   % The returned density includes the air-filled pore fraction so
   % melt-weathered ice remains on the ice path until a snow layer exists.

   %#codegen

   persistent ro_ice ro_liq ro_air
   if isempty(ro_ice)
      [ro_ice, ro_liq, ro_air] = icemodel.physicalConstant( ...
         'ro_ice', 'ro_liq', 'ro_air');
   end

   f_air_sfc = max(0, 1 - f_ice_sfc - f_liq_sfc);
   ro_sfc = ro_ice * f_ice_sfc + ro_liq * f_liq_sfc + ro_air * f_air_sfc;
end
