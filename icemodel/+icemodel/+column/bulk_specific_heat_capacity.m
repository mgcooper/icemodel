function cp = bulk_specific_heat_capacity(f_ice, f_liq, ro_bulk)
   %BULK_SPECIFIC_HEAT_CAPACITY bulk cp of ice, liquid, and air mixture.
   %
   %   cp = bulk_specific_heat_capacity(f_ice, f_liq)
   %   cp = bulk_specific_heat_capacity(f_ice, f_liq, ro_bulk)
   %
   %#codegen

   persistent cv_ice cv_liq
   if isempty(cv_ice)
      [cv_ice, cv_liq] = icemodel.physicalConstant('cv_ice', 'cv_liq');
   end

   if nargin < 3
      ro_bulk = bulk_density(f_ice, f_liq);
   end

   cp = (cv_ice * f_ice + cv_liq * f_liq) ./ ro_bulk;
end
