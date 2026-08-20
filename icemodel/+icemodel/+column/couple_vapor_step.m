function [f_ice, f_liq, budget] = couple_vapor_step( ...
      T, f_ice, f_liq, U_vap, L_vap, dz, dt, f_ice_min, f_res_por, budget)
   %COUPLE_VAPOR_STEP Apply accepted interior vapor transport once.
   %
   %  [f_ice, f_liq, budget] = icemodel.column.couple_vapor_step( ...
   %     T, f_ice, f_liq, U_vap, L_vap, dz, dt, f_ice_min, f_res_por, budget)
   %
   % Call this once per accepted substep after the surface budget closes.
   % At this point, d_liq records changes between the checkpoint and the
   % surface mass-balance call, so earlier transport would be read as melt or
   % refreezing. The vapor storage baseline covers only surface exchange, so
   % earlier transport would be scored as surface exchange.
   %
   % U_VAP is the mass flux accepted by the enthalpy solve [kg m-2 s-1],
   % positive downward, with both boundary faces closed. Its divergence
   % supplies the signed liquid-water-equivalent node increments. L_VAP is
   % the face donor latent heat from the same accepted sweep [J kg-1]: Lv
   % where the donor cell was wet at the solve state, Ls otherwise. The
   % phase route derives from it exactly (L_VAP == Lv), so the mass lands
   % in the phase whose latent heat the solve's face energy carried, with
   % no driver-side solve-state snapshot. Storage limits clamp against the
   % current fractions, because a cell can only give what it holds now.
   %
   % Each interior face is split into donor-phase removal and matching
   % receiver addition so cross-phase faces retain the donor latent heat.
   % Face 1 is closed in this interior call; the surface exchange is the
   % surface budget's business.
   %
   % Inputs
   %   T              - Column temperature at the accepted substep [K].
   %   f_ice, f_liq   - Column phase fractions after surface exchange [-].
   %   U_vap          - Accepted interior vapor mass flux [kg m-2 s-1]
   %                    (JJ+1 x 1, boundary faces zero).
   %   L_vap          - Face donor latent heat from the accepted sweep
   %                    [J kg-1] (JJ+1 x 1).
   %   dz             - Control-volume thickness [m].
   %   dt             - Substep length [s].
   %   f_ice_min      - Minimum retained ice fraction [-].
   %   f_res_por      - Residual liquid-water fraction per pore volume [-].
   %   budget         - Forcing-step budget with the post-exchange baselines
   %                    in budget.substep.
   %
   % Outputs
   %   f_ice, f_liq   - Phase fractions after the interior exchange [-].
   %   budget         - Budget with the redistribution increments.
   %
   % See also: icemodel.column.apply_vapor_transfer,
   %  icemodel.column.accumulate_redistribution_budget,
   %  icemodel.vapor.latent_enthalpy_switch
   %
   %#codegen

   persistent ro_liq Lv
   if isempty(ro_liq)
      [ro_liq, Lv] = icemodel.physicalConstant('ro_liq', 'Lv');
   end

   JJ = numel(f_ice);

   % Preserve the donor phase on every face. A node-wise net divergence cannot
   % do this across a wet/dry face because incoming and outgoing donors can
   % differ. The donor latent heat from the accepted sweep names the phase
   % exactly: both values are assigned constants, so the equality is exact.
   % The residual floor for the current state limits what a wet cell gives.
   donor_wet = L_vap(2:JJ) == Lv;
   [~, f_res] = icemodel.column.vapor_exchange_is_wet( ...
      f_ice, f_liq, f_res_por);
   d_vap_liq_nodes = zeros(JJ, 1);
   d_vap_ice_nodes = zeros(JJ, 1);
   d_vap_face = U_vap(2:JJ) * dt / ro_liq;
   d_liq_north = -d_vap_face ./ dz(1:JJ-1) .* donor_wet;
   d_liq_south = d_vap_face ./ dz(2:JJ) .* donor_wet;
   d_ice_north = -d_vap_face ./ dz(1:JJ-1) .* ~donor_wet;
   d_ice_south = d_vap_face ./ dz(2:JJ) .* ~donor_wet;
   d_vap_liq_nodes(1:JJ-1) = d_vap_liq_nodes(1:JJ-1) + d_liq_north;
   d_vap_liq_nodes(2:JJ) = d_vap_liq_nodes(2:JJ) + d_liq_south;
   d_vap_ice_nodes(1:JJ-1) = d_vap_ice_nodes(1:JJ-1) + d_ice_north;
   d_vap_ice_nodes(2:JJ) = d_vap_ice_nodes(2:JJ) + d_ice_south;

   [f_ice, f_liq] = icemodel.column.apply_vapor_transfer(f_ice, f_liq, ...
      d_vap_liq_nodes, d_vap_ice_nodes, f_ice_min, f_res);

   % Record the per-phase storage the transport moved. The baseline is the
   % post-exchange storage the surface budget wrote to budget.substep.
   budget = icemodel.column.accumulate_redistribution_budget( ...
      budget, T, f_ice, f_liq, dz);
end
