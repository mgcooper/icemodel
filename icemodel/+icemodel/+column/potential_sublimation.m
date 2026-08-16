function d_psbl = potential_sublimation(d_pevp)
   %POTENTIAL_SUBLIMATION Convert a vapor tendency to an ice fraction.
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
   % icemodel.surface.apply_surface_vapor_exchange spends the demand on
   % ice. icemodel.column.merge_thin_layers predicts whether that spending
   % takes a layer below the retained ice floor. A predictor that used a
   % different factor than the applier would flag merges the applier never
   % causes. icemodel.column.apply_vapor_transport converts a wet cell's
   % unapplied share to the same ice basis, in both of its wet branches.
   %
   % Sign convention follows d_pevp: negative is sublimation, positive is
   % deposition.
   %
   % See also: icemodel.surface.apply_surface_vapor_exchange,
   %  icemodel.column.merge_thin_layers,
   %  icemodel.surface.potential_surface_vapor_tendency
   %
   %#codegen

   persistent Ls Lv ro_ice ro_liq
   if isempty(Ls)
      [Ls, Lv, ro_ice, ro_liq] = icemodel.physicalConstant( ...
         'Ls', 'Lv', 'ro_ice', 'ro_liq');
   end

   d_psbl = d_pevp * (Lv * ro_liq) / (Ls * ro_ice);
end
