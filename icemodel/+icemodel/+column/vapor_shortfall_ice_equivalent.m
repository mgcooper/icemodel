function d_sbl_err = vapor_shortfall_ice_equivalent( ...
      d_vap_liq_unapplied, d_vap_ice_unapplied)
   %VAPOR_SHORTFALL_ICE_EQUIVALENT Convert phase shortfall for the ledger.
   %
   %  d_sbl_err = icemodel.column.vapor_shortfall_ice_equivalent( ...
   %     d_vap_liq_unapplied, d_vap_ice_unapplied)
   %
   % Both inputs are signed liquid-water-equivalent volume fractions.
   % D_SBL_ERR is the ice-fraction increment with the same latent energy:
   %
   %   ro_ice * Ls * d_sbl_err = ro_liq * (
   %      Lv * d_vap_liq_unapplied + Ls * d_vap_ice_unapplied).
   %
   % Surface and interior ledgers use this shared conversion before they
   % scale the result by cell thickness.
   %
   % See also: icemodel.column.potential_sublimation,
   %  icemodel.column.apply_vapor_transfer,
   %  icemodel.column.accumulate_redistribution_budget
   %
   %#codegen

   persistent ro_ice ro_liq
   if isempty(ro_ice)
      [ro_ice, ro_liq] = ...
         icemodel.physicalConstant('ro_ice', 'ro_liq');
   end

   d_sbl_err = icemodel.column.potential_sublimation( ...
      d_vap_liq_unapplied) ...
      + ro_liq / ro_ice * d_vap_ice_unapplied;
end
