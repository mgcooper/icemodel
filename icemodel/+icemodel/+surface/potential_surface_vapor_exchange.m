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
   % Call this once per accepted substep, BEFORE the exchange is applied. The
   % applier routes the demand on the pre-exchange state. A merge can also
   % replace the top cell outright. Reading the phase afterwards would
   % therefore classify the wrong cell.
   %
   % Accumulating the corrected fraction rather than the raw tendency is what
   % keeps a multi-substep forcing step right. The correction is linear in
   % d_pevp. Summing the corrected increments therefore equals one conversion
   % of the total whenever the phase holds, and differs correctly when it does
   % not. A step that refines through a melt or freeze transition is when the
   % substepper takes several substeps and the phase changes within one.
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
