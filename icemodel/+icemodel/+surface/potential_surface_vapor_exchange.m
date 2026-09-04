function [d_vap_liq, d_vap_ice_lwe, f_res] = potential_surface_vapor_exchange( ...
      d_pevp, f_ice, f_liq, f_res_por)
   %POTENTIAL_SURFACE_VAPOR_EXCHANGE Partition surface vapor demand.
   %
   %  [d_vap_liq, d_vap_ice_lwe, f_res] = ...
   %     icemodel.surface.potential_surface_vapor_exchange( ...
   %     d_pevp, f_ice, f_liq, f_res_por)
   %
   % D_PEVP is an energy demand expressed as a liquid-fraction increment. This
   % function partitions it into liquid and ice mass increments, both expressed
   % as liquid-fractions:
   %
   %   Lv * d_pevp = Lv * d_vap_liq + Ls * d_vap_ice_lwe.
   %
   % Evaporation first uses f_liq above F_RES, then converts any remaining
   % demand to a liquid-water-equivalent sublimation increment.
   % APPLY_VAPOR_TRANSFER limits
   % condensation and deposition by control-volume water capacity and limits
   % sublimation by F_ICE_MIN. Call this before mutating the phase fractions so
   % the partition uses the input state.
   %
   % Inputs
   %   d_pevp    - Potential surface vapor demand [-], on the Lv basis.
   %   f_ice     - Ice fraction of the top cell [-], before exchange.
   %   f_liq     - Liquid fraction of the top cell [-], before exchange.
   %   f_res_por - Residual liquid-water fraction per pore volume [-].
   %
   % Outputs
   %   d_vap_liq - Liquid-phase increment on the liquid-water basis [-].
   %   d_vap_ice_lwe - Ice-phase increment on the liquid-water basis [-].
   %   f_res     - Residual liquid fraction used for the phase decision [-].
   %
   % See also: icemodel.surface.surface_vapor_mass_flux,
   %  icemodel.column.vapor_exchange_is_wet,
   %  icemodel.surface.potential_surface_vapor_demand
   %
   %#codegen

   persistent Ls Lv
   if isempty(Ls)
      [Ls, Lv] = icemodel.physicalConstant('Ls', 'Lv');
   end

   % Identify wet cells and the residual water floor of each cell.
   [wet, f_res] = icemodel.column.vapor_exchange_is_wet( ...
      f_ice, f_liq, f_res_por);

   % Initial/default values.
   d_vap_liq = 0;
   d_vap_ice_lwe = 0;

   if wet && d_pevp < 0
      % Limit evaporation to the available liquid using Lv. Send the
      % remaining demand to ice at Ls, as a liquid fraction.
      d_vap_liq = max(d_pevp, -(f_liq - f_res));
      d_vap_ice_lwe = (d_pevp - d_vap_liq) * Lv / Ls;
   elseif wet
      d_vap_liq = d_pevp;
   else
      d_vap_ice_lwe = d_pevp * Lv / Ls;
   end
end
