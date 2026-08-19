function [d_vap_liq, d_vap_ice, f_res] = potential_surface_vapor_exchange( ...
      d_pevp, f_ice, f_liq, f_res_por)
   %POTENTIAL_SURFACE_VAPOR_EXCHANGE Partition a surface vapor demand.
   %
   %  [d_vap_liq, d_vap_ice, f_res] = ...
   %     icemodel.surface.potential_surface_vapor_exchange( ...
   %     d_pevp, f_ice, f_liq, f_res_por)
   %
   % D_PEVP is an energy demand expressed as a liquid-water-equivalent
   % fraction. This function partitions it into liquid and ice mass
   % increments, both still expressed as liquid-water-equivalent fractions.
   % Before storage limits, the partition closes the latent-energy identity
   %
   %   Lv * d_pevp = Lv * d_vap_liq + Ls * d_vap_ice.
   %
   % Wet condensation goes to liquid. Wet evaporation first exhausts liquid
   % above F_RES, then converts the remaining energy demand to ice at Lv/Ls.
   % A dry surface exchanges only with ice. Call this before mutating the
   % surface state so its phase decision matches the state being applied.
   %
   % Inputs
   %   d_pevp    - Potential surface vapor demand [-], on the Lv basis.
   %   f_ice     - Ice fraction of the top cell [-], before the exchange.
   %   f_liq     - Liquid fraction of the top cell [-], before the exchange.
   %   f_res_por - Residual liquid-water fraction per pore volume [-].
   %
   % Outputs
   %   d_vap_liq - Liquid-phase increment on the liquid-water basis [-].
   %   d_vap_ice - Ice-phase increment on the liquid-water basis [-].
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

   [wet, f_res] = icemodel.column.vapor_exchange_is_wet( ...
      f_ice, f_liq, f_res_por);

   d_vap_liq = 0;
   d_vap_ice = 0;

   if wet && d_pevp < 0
      % Spend only the available liquid at Lv. Ice receives the remaining
      % energy demand at Ls, still reported on the liquid-water mass basis.
      d_vap_liq = max(d_pevp, -(f_liq - f_res));
      d_vap_ice = (d_pevp - d_vap_liq) * Lv / Ls;
   elseif wet
      d_vap_liq = d_pevp;
   else
      d_vap_ice = d_pevp * Lv / Ls;
   end
end
