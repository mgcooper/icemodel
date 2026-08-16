function d_liq_max = max_liquid_fraction_change(f_ice, f_liq)
   %MAX_LIQUID_FRACTION_CHANGE Largest f_liq increase a control volume takes.
   %
   %  d_liq_max = icemodel.column.max_liquid_fraction_change(f_ice, f_liq)
   %
   % Returns the largest increase in f_liq a control volume accepts before its
   % total water reaches f_wat_max:
   %
   %   d_liq_max = ro_ice / ro_liq * (1 - f_ice) - f_liq
   %             = f_wat_max - f_wat
   %
   % on the icemodel.column.water_fraction basis, where
   % f_wat = f_liq + f_ice * ro_ice / ro_liq and f_wat_max = ro_ice / ro_liq.
   % Both sides are liquid-water-equivalent volume fractions.
   %
   % The pore volume 1 - f_ice is scaled to water equivalent as if it were
   % ice. Write P = 1 - f_ice and r = ro_ice / ro_liq. The bound is r * P,
   % which falls short of P by (1 - r) * P. That shortfall is the room the
   % liquid needs to expand if it refreezes. It keeps
   % f_ice + f_liq * ro_liq / ro_ice from exceeding 1.
   %
   % Read 1 - f_ice as open porosity in snow. In bubbly glacier ice it is the
   % total non-ice volume. Some of that volume is closed to liquid there, so
   % d_liq_max is an upper bound. The model carries no bubble fraction, so
   % f_ice and f_liq are the only phase state a cell can compute this from.
   %
   % One function owns this limit, so the surface exchange and the interior
   % transport cannot cap condensation differently.
   %
   % Inputs
   %   f_ice    - Ice fraction [-], scalar or array.
   %   f_liq    - Liquid-water fraction [-], same shape as f_ice.
   %
   % Output
   %   d_liq_max - Largest f_liq increase the cell accepts [-]. Negative where
   %               a cell is already over f_wat_max; callers floor or clamp.
   %
   % See also: icemodel.surface.apply_surface_vapor_exchange,
   %  icemodel.column.apply_vapor_transport,
   %  icemodel.column.infiltration,
   %  icemodel.column.assert_max_water
   %
   %#codegen

   persistent ro_ice ro_liq
   if isempty(ro_ice)
      [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   end

   d_liq_max = ro_ice ./ ro_liq .* (1.0 - f_ice) - f_liq;
end
