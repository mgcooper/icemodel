function [f_ice, f_liq, d_vap_liq_unapplied, d_vap_ice_unapplied] = ...
      apply_vapor_transfer(f_ice, f_liq, d_vap_liq_nodes, ...
      d_vap_ice_nodes, f_ice_min, f_res)
   %APPLY_VAPOR_TRANSFER Apply phase-resolved vapor increments to the column.
   %
   %  [f_ice, f_liq, d_vap_liq_unapplied, d_vap_ice_unapplied] = ...
   %     icemodel.column.apply_vapor_transfer(f_ice, f_liq, ...
   %     d_vap_liq_nodes, d_vap_ice_nodes, f_ice_min, f_res)
   %
   % Both input increments and both unapplied outputs are signed
   % liquid-water-equivalent volume fractions. Positive increments add mass;
   % negative increments remove it. Phase selection and energy-to-mass
   % conversion happen upstream. This function only mutates the selected
   % phase under the shared storage limits: liquid condensation is capped by
   % pore capacity, liquid evaporation stops at F_RES, ice deposition is
   % capped by representable water-storage capacity, and ice sublimation
   % stops at F_ICE_MIN.
   %
   % See also: icemodel.surface.apply_surface_vapor_exchange,
   %  icemodel.column.couple_vapor_step
   %
   %#codegen

   persistent ro_ice ro_liq
   if isempty(ro_ice)
      [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   end

   d_vap_liq_unapplied = zeros(size(d_vap_liq_nodes));
   d_vap_ice_unapplied = zeros(size(d_vap_ice_nodes));

   % Apply all liquid increments first. This preserves the scalar loop's
   % liquid-before-ice capacity ordering while avoiding one function call per
   % column node on every accepted substep.
   positive_liq = d_vap_liq_nodes > 0;
   negative_liq = d_vap_liq_nodes < 0;
   capacity_liq = max(icemodel.column.max_liquid_fraction_change( ...
      f_ice, f_liq), 0);
   applied_liq = zeros(size(d_vap_liq_nodes));
   applied_liq(positive_liq) = min( ...
      d_vap_liq_nodes(positive_liq), capacity_liq(positive_liq));
   available_liq = max(f_liq - f_res, 0);
   applied_liq(negative_liq) = max( ...
      d_vap_liq_nodes(negative_liq), -available_liq(negative_liq));
   d_vap_liq_unapplied(positive_liq) = ...
      d_vap_liq_nodes(positive_liq) - applied_liq(positive_liq);
   d_vap_liq_unapplied(negative_liq) = ...
      d_vap_liq_nodes(negative_liq) - applied_liq(negative_liq);
   f_liq = f_liq + applied_liq;

   % Apply ice increments after liquid mutation. Positive ice transfer stays
   % on the LWE basis so residual liquid retains its expansion room; negative
   % ice transfer uses the stored ice-volume basis and converts its shortfall
   % back to LWE.
   positive_ice = d_vap_ice_nodes > 0;
   negative_ice = d_vap_ice_nodes < 0;
   capacity_ice = max(icemodel.column.max_liquid_fraction_change( ...
      f_ice, f_liq), 0);
   applied_ice_lwe = zeros(size(d_vap_ice_nodes));
   applied_ice_lwe(positive_ice) = min( ...
      d_vap_ice_nodes(positive_ice), capacity_ice(positive_ice));
   d_vap_ice_unapplied(positive_ice) = ...
      d_vap_ice_nodes(positive_ice) - applied_ice_lwe(positive_ice);
   f_ice = f_ice + applied_ice_lwe * ro_liq / ro_ice;

   d_ice = d_vap_ice_nodes * ro_liq / ro_ice;
   available_ice = max(f_ice - f_ice_min, 0);
   applied_ice = zeros(size(d_vap_ice_nodes));
   applied_ice(negative_ice) = max( ...
      d_ice(negative_ice), -available_ice(negative_ice));
   d_vap_ice_unapplied(negative_ice) = ...
      (d_ice(negative_ice) - applied_ice(negative_ice)) * ro_ice / ro_liq;
   f_ice = f_ice + applied_ice;
end
