function [f_ice, f_liq, d_vap_liq_unapplied, d_vap_ice_unapplied_lwe] = ...
      apply_vapor_transfer(f_ice, f_liq, d_vap_liq_dmd, d_vap_ice_dmd_lwe, ...
      f_ice_min, f_res)
   %APPLY_VAPOR_TRANSFER Apply vapor phase change increments to the column.
   %
   %  [f_ice, f_liq, d_vap_liq_unapplied, d_vap_ice_unapplied_lwe] = ...
   %     icemodel.column.apply_vapor_transfer(f_ice, f_liq, ...
   %     d_vap_liq_dmd, d_vap_ice_dmd_lwe, f_ice_min, f_res)
   %
   % Converts potential vapor phase change increments to applied increments,
   % updates the ice and liquid water fractions, and returns unapplied
   % demand. The function is called by the surface vapor exchange and
   % subsurface vapor transport models.
   %
   % The input and unapplied demand are signed liquid-water-equivalent volume
   % fractions. Positive values add water; negative values remove it. Liquid
   % condensation stops at control-volume water capacity, evaporation stops at
   % F_RES, ice deposition stops at control-volume water capacity, and
   % sublimation stops at F_ICE_MIN.
   %
   % See also: icemodel.surface.apply_surface_vapor_exchange,
   %  icemodel.column.couple_vapor_step
   %
   %#codegen

   persistent ro_ice ro_liq
   if isempty(ro_ice)
      [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   end

   % Initialize vectors. Start with no unapplied vapor transfer.
   applied_liq = zeros(size(d_vap_liq_dmd));
   applied_ice = zeros(size(d_vap_ice_dmd_lwe));
   applied_ice_lwe = zeros(size(d_vap_ice_dmd_lwe));
   d_vap_liq_unapplied = zeros(size(d_vap_liq_dmd));
   d_vap_ice_unapplied_lwe = zeros(size(d_vap_ice_dmd_lwe));

   % Transfer liquid first then ice, same as apply_surface_vapor_exchange.

   % Addition and removal use different storage limits, so find cells with
   % positive and negative demand to evaluate the limits separately.
   positive_liq = d_vap_liq_dmd > 0;
   negative_liq = d_vap_liq_dmd < 0;

   % Limit condensation by the control-volume water capacity.
   capacity_liq = max(icemodel.column.max_liquid_fraction_change( ...
      f_ice, f_liq), 0);
   applied_liq(positive_liq) = min( ...
      d_vap_liq_dmd(positive_liq), capacity_liq(positive_liq));

   % Limit evaporation to water above the residual liquid f_res.
   available_liq = max(f_liq - f_res, 0);
   applied_liq(negative_liq) = max( ...
      d_vap_liq_dmd(negative_liq), -available_liq(negative_liq));

   % Track the unapplied liquid change.
   d_vap_liq_unapplied(positive_liq) = ...
      d_vap_liq_dmd(positive_liq) - applied_liq(positive_liq);
   d_vap_liq_unapplied(negative_liq) = ...
      d_vap_liq_dmd(negative_liq) - applied_liq(negative_liq);

   % Update the liquid fraction.
   f_liq = f_liq + applied_liq;

   % Transfer ice after liquid.
   positive_ice = d_vap_ice_dmd_lwe > 0;
   negative_ice = d_vap_ice_dmd_lwe < 0;

   % Limit deposition by the control-volume water capacity.
   capacity_ice = max(icemodel.column.max_liquid_fraction_change( ...
      f_ice, f_liq), 0);
   applied_ice_lwe(positive_ice) = min( ...
      d_vap_ice_dmd_lwe(positive_ice), capacity_ice(positive_ice));

   % Add deposition to f_ice in ice-fraction units and convert sublimation
   % demand to ice-fraction units before applying it.
   f_ice = f_ice + applied_ice_lwe * ro_liq / ro_ice;
   d_ice = d_vap_ice_dmd_lwe * ro_liq / ro_ice;

   % Limit sublimation to ice above f_ice_min.
   available_ice = max(f_ice - f_ice_min, 0);
   applied_ice(negative_ice) = max( ...
      d_ice(negative_ice), -available_ice(negative_ice));

   % Track unapplied deposition and sublimation, converting sublimation back to
   % liquid-fraction units.
   d_vap_ice_unapplied_lwe(positive_ice) = ...
      d_vap_ice_dmd_lwe(positive_ice) - applied_ice_lwe(positive_ice);
   d_vap_ice_unapplied_lwe(negative_ice) = ...
      (d_ice(negative_ice) - applied_ice(negative_ice)) * ro_ice / ro_liq;

   % Update the ice fraction.
   f_ice = f_ice + applied_ice;
end
