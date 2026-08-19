function d_psbl = potential_sublimation(d_pevp)
   %POTENTIAL_SUBLIMATION Convert a vapor demand to an ice fraction.
   %
   %  d_psbl = icemodel.column.potential_sublimation(d_pevp)
   %
   % D_PEVP is the potential vapor-driven change expressed as a liquid-water
   % volume fraction. D_PSBL is the same latent-heat demand expressed as an
   % ice volume fraction. The surface energy balance produces one latent heat
   % flux Qe, and the two fractions are the two ways to spend it:
   %
   %   d_pevp = Qe / (Lv * ro_liq) * dt / dz
   %   d_psbl = Qe / (Ls * ro_ice) * dt / dz
   %
   % so d_psbl = d_pevp * (Lv * ro_liq) / (Ls * ro_ice).
   %
   % Three callers need this conversion and must agree by construction.
   % icemodel.surface.apply_surface_vapor_exchange converts phase-resolved
   % surface shortfall to the legacy ice-fraction record.
   % icemodel.column.merge_thin_layers predicts whether surface spending takes
   % a layer below the retained ice floor. icemodel.column.couple_vapor_step
   % converts liquid-phase interior shortfall to the same ice-fraction energy
   % basis for the redistribution ledger.
   %
   % Sign convention follows d_pevp: negative is sublimation, positive is
   % deposition.
   %
   % See also: icemodel.surface.apply_surface_vapor_exchange,
   %  icemodel.column.merge_thin_layers,
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
