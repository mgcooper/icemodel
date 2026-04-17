function ro = bulk_density(f_ice, f_liq)
   %BULK_DENSITY bulk density of ice, liquid, and air mixture
   %
   %   ro = bulk_density(f_ice, f_liq)
   %
   %#codegen

   persistent ro_ice ro_liq ro_air
   if isempty(ro_ice)
      [ro_ice, ro_liq, ro_air] = icemodel.physicalConstant( ...
         'ro_ice', 'ro_liq', 'ro_air');
   end

   ro = ro_ice * f_ice + ro_liq * f_liq + ro_air * (1.0 - f_liq - f_ice);
end
