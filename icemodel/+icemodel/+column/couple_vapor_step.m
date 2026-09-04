function [f_ice, f_liq, d_vap_liq, d_vap_ice, budget] = ...
      couple_vapor_step(f_ice, f_liq, U_vap, L_vap, d_vap_liq, ...
      d_vap_ice, dz, dt, f_ice_min, f_res_por, budget)
   %COUPLE_VAPOR_STEP Apply subsurface vapor transport across cells.
   %
   %  [f_ice, f_liq, d_vap_liq, d_vap_ice, budget] = ...
   %     icemodel.column.couple_vapor_step(f_ice, f_liq, U_vap, L_vap, ...
   %     d_vap_liq, d_vap_ice, dz, dt, f_ice_min, f_res_por, budget)
   %
   % This function applies subsurface vapor transport driven by the
   % temperature-dependent subsurface saturated vapor pressure gradient at the
   % end of each substep and tracks the accumulated applied phase changes. The
   % enthalpy solve calls icemodel.column.vapor_transport_terms and returns its
   % U_vap and L_vap outputs through the coupler. The terms applied here are
   % therefore from the accepted solve. Vapor transport changes the phase
   % fractions, so the driver calls budget_surface_mass_balance to record
   % melt/freeze phase change before calling this function.
   %
   % U_VAP is the mass flux from the enthalpy solve [kg m-2 s-1], positive
   % downward, with both boundary faces closed. Positive flux removes mass from
   % the cell above and adds it below and negative flux vise versa.
   % Dividing each change by cell thickness gives the signed
   % liquid-water-equivalent increments. L_VAP is the donor-cell latent heat
   % from the solve [J kg-1]: Lv for a wet donor and Ls otherwise. The latent
   % heat selects whether liquid or ice transfers, so the mass is in the same
   % phase as the latent heat used for the solve's face energy. Control-volume
   % water capacity, residual liquid water, and the minimum ice fraction limit
   % the applied transfer. The function subtracts rejected increments from the
   % demand and records the applied storage change.
   %
   % Inputs
   %   f_ice, f_liq   - Column phase fractions after surface exchange [-].
   %   U_vap          - Interior vapor mass flux [kg m-2 s-1]
   %                    (JJ+1 x 1, boundary faces zero).
   %   L_vap          - Face donor latent heat from the solve
   %                    [J kg-1] (JJ+1 x 1).
   %   d_vap_liq      - Accumulated vapor-driven liquid increments per
   %                    cell [-] (exchange plus transport; written out as
   %                    ice2.df_vap_liq).
   %   d_vap_ice      - Accumulated vapor-driven ice increments per
   %                    cell [-], ice fraction basis (written out as
   %                    ice2.df_vap_ice).
   %   dz             - Control-volume thickness [m].
   %   dt             - Substep length [s].
   %   f_ice_min      - Minimum allowed ice fraction [-].
   %   f_res_por      - Residual liquid-water fraction per pore volume [-].
   %   budget         - Forcing-step budget
   %                    (icemodel.column.initialize_budget_state).
   %
   % Outputs
   %   f_ice, f_liq   - Phase fractions after the interior transport [-].
   %   d_vap_liq      - Accumulator with this substep's applied transport
   %                    added.
   %   d_vap_ice      - Accumulator with this substep's applied transport
   %                    added, ice fraction units.
   %   budget         - Budget with the transport increments added.
   %
   % See also: icemodel.column.apply_vapor_transfer,
   %  icemodel.column.accumulate_vapor_transport,
   %  icemodel.column.vapor_transport_terms
   %
   %#codegen

   persistent ro_ice ro_liq Lv
   if isempty(ro_liq)
      [ro_ice, ro_liq, Lv] = ...
         icemodel.physicalConstant('ro_ice', 'ro_liq', 'Lv');
   end

   JJ = numel(f_ice);

   % Initialize vectors for the demand increments (_dmd = demand).
   d_vap_liq_dmd = zeros(JJ, 1);
   d_vap_ice_dmd_lwe = zeros(JJ, 1);

   % The residual-water fraction floor limits what a wet cell can give.
   f_res = icemodel.column.residual_water_fraction(f_ice, f_liq, f_res_por);

   % Use the donor latent heat on each face.
   donor_wet = L_vap(2:JJ) == Lv;

   % Convert face flux to liquid-water-equivalent depth. Remove it from the
   % north cell and add it to the south cell using the donor phase.
   d_vap_faces = U_vap(2:JJ) * dt / ro_liq;
   d_liq_north = -d_vap_faces ./ dz(1:JJ-1) .* donor_wet;
   d_liq_south = d_vap_faces ./ dz(2:JJ) .* donor_wet;
   d_ice_north = -d_vap_faces ./ dz(1:JJ-1) .* ~donor_wet;
   d_ice_south = d_vap_faces ./ dz(2:JJ) .* ~donor_wet;

   % Add the north- and south-face contributions in each cell.
   d_vap_liq_dmd(1:JJ-1) = d_vap_liq_dmd(1:JJ-1) + d_liq_north;
   d_vap_liq_dmd(2:JJ) = d_vap_liq_dmd(2:JJ) + d_liq_south;
   d_vap_ice_dmd_lwe(1:JJ-1) = d_vap_ice_dmd_lwe(1:JJ-1) + d_ice_north;
   d_vap_ice_dmd_lwe(2:JJ) = d_vap_ice_dmd_lwe(2:JJ) + d_ice_south;

   % Apply the partitioned solid/liquid demand to f_ice/f_liq and return any
   % unapplied demand rejected by water capacity, residual liquid f_res, or
   % minimum ice fraction f_ice_min.
   [f_ice, f_liq, d_vap_liq_unapplied, d_vap_ice_unapplied_lwe] = ...
      icemodel.column.apply_vapor_transfer(f_ice, f_liq, ...
      d_vap_liq_dmd, d_vap_ice_dmd_lwe, f_ice_min, f_res);

   % Subtract the unapplied amounts from the demand to get the applied changes.
   d_vap_liq_applied = d_vap_liq_dmd - d_vap_liq_unapplied;
   d_vap_ice_applied_lwe = d_vap_ice_dmd_lwe - d_vap_ice_unapplied_lwe;

   % Convert the ice change to ice-fraction units for d_vap_ice and the budget.
   d_vap_ice_applied = d_vap_ice_applied_lwe * ro_liq / ro_ice;

   % Budget interior vapor transport. d_vap_liq and d_vap_ice accumulate surface
   % vapor exchange in budget_surface_mass_balance and interior transport here.
   d_vap_liq = d_vap_liq + d_vap_liq_applied;
   d_vap_ice = d_vap_ice + d_vap_ice_applied;

   % Add the applied transport to the forcing-step budget.
   budget = icemodel.column.accumulate_vapor_transport( ...
      budget, d_vap_liq_applied, d_vap_ice_applied, dz);
end
