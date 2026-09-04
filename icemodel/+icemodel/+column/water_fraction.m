function f_wat = water_fraction(f_ice, f_liq)
   %WATER_FRACTION Compute the total volumetric water fraction.
   %
   %  f_wat = icemodel.column.water_fraction(f_ice, f_liq)
   %
   % Returns the total water fraction for a control volume, expressed as the
   % sum of liquid-water volume plus the liquid-water-equivalent volume of the
   % ice fraction:
   %
   %   f_wat = f_liq + f_ice * ro_ice / ro_liq
   %
   % See also: icemodel.column.bulk_enthalpy,
   %  icemodel.column.max_liquid_fraction_change
   %
   %#codegen

   persistent ro_ice ro_liq
   if isempty(ro_ice)
      [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   end

   f_wat = f_liq + f_ice * ro_ice / ro_liq;
end
