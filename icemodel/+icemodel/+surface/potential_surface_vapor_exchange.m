function d_vap_sfc = potential_surface_vapor_exchange( ...
      d_pevp, f_ice, f_liq, f_res_por)
   %POTENTIAL_SURFACE_VAPOR_EXCHANGE Phase-correct a surface vapor demand.
   %
   %  d_vap_sfc = icemodel.surface.potential_surface_vapor_exchange( ...
   %     d_pevp, f_ice, f_liq, f_res_por)
   %
   % D_PEVP is an ENERGY demand wearing the units of a liquid-water volume
   % fraction: the surface energy balance forms it as Qe / (Lv * ro_liq) *
   % dt / dz whatever phase the surface is in. The mass that demand moves
   % depends on the phase supplying it, because the same energy sublimates
   % less mass than it evaporates. This returns the demand rescaled to the
   % liquid-water volume fraction of the mass it actually moves:
   %
   %   d_vap_sfc = d_pevp * Lv / L_sfc
   %
   % A wet surface exchanges at Lv, so the latent heats cancel and the demand
   % passes through. A dry surface exchanges at Ls, so the fraction shrinks by
   % Lv / Ls.
   %
   % Call this with the pre-exchange state: the applier routes the demand
   % on that state, so a later read would classify the wrong cell.
   %
   % No production path calls this. The driver accumulates the realized
   % exchange the applier returns, which already carries the phase the
   % limits actually spent. This converter remains the single owner of the
   % demand-to-mass correction for the planned unification that converts
   % the surface demand upstream and routes it through one column applier.
   %
   % One latent heat per substep is still an approximation. A demand larger
   % than the mobile liquid spends the remainder on ice inside the same
   % applier call. This charges all of it to the phase the cell started in.
   % Bead icemodel-6nv measures that.
   %
   % Inputs
   %   d_pevp    - Potential surface vapor tendency [-], on the Lv basis.
   %   f_ice     - Ice fraction of the top cell [-], before the exchange.
   %   f_liq     - Liquid fraction of the top cell [-], before the exchange.
   %   f_res_por - Residual liquid-water fraction per pore volume [-].
   %
   % Output
   %   d_vap_sfc - Liquid-water volume fraction of the mass the demand moves.
   %
   % See also: icemodel.surface.surface_vapor_mass_flux,
   %  icemodel.column.vapor_exchange_is_wet,
   %  icemodel.surface.potential_surface_vapor_tendency
   %
   %#codegen

   persistent Ls Lv
   if isempty(Ls)
      [Ls, Lv] = icemodel.physicalConstant('Ls', 'Lv');
   end

   % One owner decides wet versus dry, so this cannot disagree with the
   % applier that spends the demand.
   if icemodel.column.vapor_exchange_is_wet(f_ice, f_liq, f_res_por)
      L_sfc = Lv;
   else
      L_sfc = Ls;
   end

   d_vap_sfc = d_pevp * Lv / L_sfc;
end
