function d_psbl = potential_sublimation(d_pevp)
   %POTENTIAL_SUBLIMATION Convert surface vapor demand to ice fraction.
   %
   %  d_psbl = icemodel.column.potential_sublimation(d_pevp)
   %
   % D_PEVP is the potential surface vapor exchange expressed as a liquid-
   % water volume fraction. D_PSBL is the same latent-energy demand expressed
   % as an ice volume fraction:
   %
   %   d_pevp = Qe / (Lv * ro_liq) * dt / dz
   %   d_psbl = Qe / (Ls * ro_ice) * dt / dz
   %
   % so d_psbl = d_pevp * (Lv * ro_liq) / (Ls * ro_ice).
   %
   % Negative values are sublimation and positive values are deposition.
   %
   % See also: icemodel.column.merge_thin_layers,
   %  icemodel.surface.potential_surface_vapor_demand
   %
   %#codegen

   persistent Ls Lv ro_ice ro_liq
   if isempty(Ls)
      [Ls, Lv, ro_ice, ro_liq] = icemodel.physicalConstant( ...
         'Ls', 'Lv', 'ro_ice', 'ro_liq');
   end

   d_psbl = d_pevp * (Lv * ro_liq) / (Ls * ro_ice);
end
