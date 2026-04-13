function tf = assert_max_water(f_ice, f_liq)

   % Load the ratio ro_ice / ro_liq ("ro_ice_water_equivalent")
   persistent ro_iwe
   if isempty(ro_iwe)
      ro_iwe = icemodel.physicalConstant('ro_iwe');
   end

   tf = all(icemodel.column.water_fraction(f_ice, f_liq) <= ro_iwe + eps);
end
