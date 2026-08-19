function [f_ice, f_liq, d_rof, d_sbl_err, d_applied] = ...
      apply_surface_vapor_exchange( ...
      f_ice, f_liq, d_rof, d_pevp, f_ice_min, f_res_por)
   %APPLY_SURFACE_VAPOR_EXCHANGE Apply surface vapor energy demand.
   %
   %  [f_ice, f_liq, d_rof, d_sbl_err, d_applied] = ...
   %     icemodel.surface.apply_surface_vapor_exchange( ...
   %     f_ice, f_liq, d_rof, d_pevp, f_ice_min, f_res_por)
   %
   % D_PEVP is the surface latent-energy demand expressed as a
   % liquid-water-equivalent top-cell fraction. The potential-exchange
   % helper partitions that demand exactly between liquid at Lv and ice at
   % Ls. The shared column transfer then applies both mass increments under
   % the phase storage limits.
   %
   % Positive liquid increment that the top cell cannot hold has already
   % crossed the surface and becomes runoff. Other unapplied demand remains
   % in D_SBL_ERR on the ice-fraction energy basis. D_APPLIED is the
   % signed liquid-water-equivalent mass that crossed the surface, including
   % condensation overflow but excluding rejected deposition or unsatisfied
   % evaporation and sublimation.
   %
   % Inputs
   %   f_ice, f_liq - Column phase fractions [-].
   %   d_rof        - Running condensation runoff fraction [-].
   %   d_pevp       - Surface vapor energy demand on the Lv basis [-].
   %   f_ice_min    - Minimum retained ice fraction [-].
   %   f_res_por    - Residual liquid fraction per pore volume [-].
   %
   % Outputs
   %   f_ice, f_liq - Updated column phase fractions [-].
   %   d_rof        - Running condensation runoff fraction [-].
   %   d_sbl_err    - Unapplied demand in ice-fraction energy units [-].
   %   d_applied    - Realized exchange on the liquid-water basis [-].
   %
   % See also: icemodel.surface.potential_surface_vapor_exchange,
   %  icemodel.column.apply_vapor_transfer
   %
   %#codegen

   persistent ro_ice ro_liq
   if isempty(ro_ice)
      [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   end

   % Retain the top-cell state so realized mass comes from the mutation the
   % limits actually produced.
   f_ice_top = f_ice(1);
   f_liq_top = f_liq(1);

   % Partition the energy demand before applying either phase increment.
   [d_vap_liq, d_vap_ice, f_res] = ...
      icemodel.surface.potential_surface_vapor_exchange( ...
      d_pevp, f_ice_top, f_liq_top, f_res_por);

   % A solved surface cell can already be below the remesh floor. In that
   % case, allow its existing ice to sublimate before merge_thin_layers runs.
   surface_f_ice_min = f_ice_min;
   if f_ice_top < f_ice_min
      surface_f_ice_min = 0;
   end

   [f_ice(1), f_liq(1), d_vap_liq_unapplied, ...
      d_vap_ice_unapplied] = icemodel.column.apply_vapor_transfer( ...
      f_ice_top, f_liq_top, d_vap_liq, d_vap_ice, ...
      surface_f_ice_min, f_res);

   % Surface liquid condensation overflow crossed the boundary and runs off.
   overflow = max(d_vap_liq_unapplied, 0);
   d_rof = d_rof + overflow;
   d_vap_liq_unapplied = min(d_vap_liq_unapplied, 0);

   % Preserve the ledger's ice-fraction energy basis.
   d_sbl_err = zeros(numel(f_ice), 1);
   d_sbl_err(1) = icemodel.column.vapor_shortfall_ice_equivalent( ...
      d_vap_liq_unapplied, d_vap_ice_unapplied);

   % Overflow crossed before runoff, so include it in realized exchange.
   d_applied = (f_liq(1) - f_liq_top) ...
      + (f_ice(1) - f_ice_top) * ro_ice / ro_liq + overflow;
end
