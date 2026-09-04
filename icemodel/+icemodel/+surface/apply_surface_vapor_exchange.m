function [f_ice, f_liq, d_rof, d_vap_liq, d_vap_ice, d_vap] = ...
      apply_surface_vapor_exchange(f_ice, f_liq, d_rof, d_pevp, dz, ...
      f_ice_min, f_res_por)
   %APPLY_SURFACE_VAPOR_EXCHANGE Apply surface vapor energy demand.
   %
   %  [f_ice, f_liq, d_rof, d_vap_liq, d_vap_ice, d_vap] = ...
   %     icemodel.surface.apply_surface_vapor_exchange( ...
   %     f_ice, f_liq, d_rof, d_pevp, dz, f_ice_min, f_res_por)
   %
   % D_PEVP is the surface latent-energy demand expressed as a
   % liquid-water-equivalent fraction of the top cell.
   % icemodel.surface.potential_surface_vapor_exchange partitions D_PEVP between
   % liquid at Lv and ice at Ls. apply_vapor_transfer applies those partitioned
   % increments limited by the liquid or ice available in the cell.
   %
   % Condensation beyond the top-cell water capacity becomes runoff. Evaporation
   % and sublimation continue into deeper cells until the demand is satisfied or
   % every cell reaches its retained liquid or ice limit. The unsatisfied demand
   % is scaled by dz(j)/dz(j+1) before it is applied to the next cell to
   % preserve areal latent energy on a nonuniform mesh.
   %
   % The function can leave part of an evaporation or sublimation demand
   % unapplied when the column has no removable water. A full top cell can also
   % reject ice deposition because the model does not add a new surface cell and
   % deposition does not cascade into deeper cells. D_VAP is the applied
   % liquid-fraction-equivalent mass that actually crossed the surface,
   % including condensation overflow.
   %
   % Outputs
   %   f_ice, f_liq    - Updated column phase fractions [-].
   %   d_rof           - Running condensation runoff fraction [-].
   %   d_vap_liq       - Applied liquid storage increments per cell [-],
   %                     liquid-water-equivalent fraction basis. Overflow
   %                     is credited to runoff (d_rof) without delay.
   %   d_vap_ice       - Applied ice increments per cell [-], ice fraction
   %                     basis (multiply by ro_ice/ro_liq for water
   %                     equivalent).
   %   d_vap           - Applied exchange on the liquid-water basis,
   %                     expressed as a fraction of the top cell [-].
   %
   % See also: icemodel.surface.potential_surface_vapor_exchange,
   %  icemodel.column.apply_vapor_transfer
   %
   %#codegen

   persistent ro_ice ro_liq Ls Lv
   if isempty(ro_ice)
      [ro_ice, ro_liq, Ls, Lv] = ...
         icemodel.physicalConstant('ro_ice', 'ro_liq', 'Ls', 'Lv');
   end

   % Initial values.
   JJ = numel(f_ice);
   d_vap_liq = zeros(JJ, 1);
   d_vap_ice = zeros(JJ, 1);
   overflow = 0.0;

   % Partition the top-cell energy demand before changing either phase.
   [d_vap_liq_dmd, d_vap_ice_dmd_lwe, f_res] = ...
      icemodel.surface.potential_surface_vapor_exchange( ...
      d_pevp, f_ice(1), f_liq(1), f_res_por);

   % Apply the demand to cell 1, then continue unsatisfied evaporation or
   % sublimation into the cells below. Deposition applies only to cell 1.
   for j = 1:JJ
      f_ice_old = f_ice(j);
      f_liq_old = f_liq(j);

      % Let a cell below F_ICE_MIN sublimate its remaining ice before remeshing.
      cell_f_ice_min = f_ice_min;
      if f_ice_old < f_ice_min
         cell_f_ice_min = 0;
      end

      % Apply the partitioned solid/liquid demand to f_ice/f_liq and get the
      % unapplied remaining demand.
      [f_ice(j), f_liq(j), d_vap_liq_unapplied, d_vap_ice_unapplied_lwe] = ...
         icemodel.column.apply_vapor_transfer(f_ice_old, f_liq_old, ...
         d_vap_liq_dmd, d_vap_ice_dmd_lwe, cell_f_ice_min, f_res);

      % Send unapplied liquid demand above the top cell's capacity to runoff.
      if j == 1
         overflow = max(d_vap_liq_unapplied, 0);
         d_rof = d_rof + overflow;
         d_vap_liq_unapplied = min(d_vap_liq_unapplied, 0);
      end

      % Difference the new and old state to get the change increments.
      d_vap_liq(j) = (f_liq(j) - f_liq_old);
      d_vap_ice(j) = (f_ice(j) - f_ice_old);

      % Break on rejected surface ice deposition (no cascade to deeper cells).
      % Positive ice demand occurs only for a dry surface, so
      % d_vap_ice_unapplied_lwe > 0 means d_vap_liq_unapplied == 0
      if d_vap_ice_unapplied_lwe > 0
         break
      end

      % Break when no unapplied demand remains or the bottom cell is reached.
      if d_vap_liq_unapplied == 0 && d_vap_ice_unapplied_lwe == 0
         break
      end
      if j == JJ
         break
      end

      % Update the remaining unsatisfied demand on the Lv energy basis and
      % scale the demand increments to the thickness of the next cell.
      d_pevp_remainder = (d_vap_liq_unapplied ...
         + d_vap_ice_unapplied_lwe * Ls / Lv) * dz(j) / dz(j + 1);
      [d_vap_liq_dmd, d_vap_ice_dmd_lwe, f_res] = ...
         icemodel.surface.potential_surface_vapor_exchange( ...
         d_pevp_remainder, f_ice(j + 1), f_liq(j + 1), f_res_por);
   end

   % Express the applied column mass as a liquid-water fraction of cell 1.
   % Include condensation runoff because it crossed the surface.
   d_vap = (sum(d_vap_liq .* dz) ...
      + sum(d_vap_ice .* dz) * ro_ice / ro_liq) / dz(1) + overflow;
end
