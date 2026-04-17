function tf = assert_max_water(f_ice, f_liq)
   %ASSERT_MAX_WATER Assert that water fraction does not exceed the maximum.
   %
   %  tf = icemodel.column.assert_max_water(f_ice, f_liq)
   %
   % Returns true if all elements satisfy f_wat <= ro_ice/ro_liq + eps.
   %
   % Physics Background
   % ------------------
   % In the current formulation the ice column is initialized to glacier ice,
   % which has closed air bubbles. The pore space is therefore treated as
   % sealed to liquid infiltration from above. When all the ice in a control
   % volume melts, the volume originally occupied by ice becomes liquid water,
   % so the water equivalent is
   %
   %   f_wat_max = ro_ice / ro_liq   (≈ 0.917)
   %
   % This is the quantity checked here. The assertion is f_wat <= f_wat_max
   % where f_wat = f_liq + f_ice * ro_ice / ro_liq (see water_fraction).
   %
   % Top-node condensation
   % ---------------------
   % The top control volume can receive condensation through
   % icemodel.surface.apply_surface_vapor_mass_change. In that path the
   % condensation capacity for a wet node is
   %
   %   d_liq_max = ro_ice/ro_liq * (1 - f_ice_top) - f_liq_top
   %
   % which is exactly the remaining headroom below f_wat_max. Condensation
   % exceeding that cap is tracked in d_rof rather than stored, so the
   % assertion should still hold for the top node even when condensation is
   % active. The check here uses `all` so a single failing node is caught.
   %
   % Melt-zone context (meltzone_transform)
   % ---------------------------------------
   % Inside meltzone_transform the liquid fraction change is bounded by
   %
   %   d_liq_max = f_wat_max - f_liq_old
   %             = ro_ice / ro_liq - f_liq(iM)
   %
   % Note that f_wat here equals ro_ice/ro_liq only because the melt-zone
   % calculation preserves total water mass: f_wat does not change as ice
   % melts, it is simply partitioned differently between f_ice and f_liq.
   % The assertion checks the result of that partitioning.
   %
   % Infiltration / snow
   % --------------------
   % When liquid infiltration from above is enabled the correct f_wat_max
   % changes because additional water can be added to the pore space beyond
   % the closed-bubble equivalent. In that regime:
   %
   %   - If f_wat is allowed to reach 1, refreezing of excess liquid can
   %     still cause f_liq to temporarily exceed ro_ice/ro_liq if the
   %     freeze path does not first cap the liquid fraction.
   %   - The correct assertion in an open-pore / snow physics context would
   %     be f_wat <= 1, not f_wat <= ro_ice/ro_liq.
   %   - Until infiltration is implemented this stricter glacier-ice bound
   %     is appropriate. If infiltration is added, revisit this check.
   %
   % See also: icemodel.column.water_fraction,
   %           icemodel.surface.apply_surface_vapor_mass_change,
   %           icemodel.column.meltzone_transform
   %
   %#codegen

   % Load the ratio ro_ice / ro_liq ("ro_ice_water_equivalent")
   persistent ro_iwe
   if isempty(ro_iwe)
      ro_iwe = icemodel.physicalConstant('ro_iwe');
   end

   tf = all(icemodel.column.water_fraction(f_ice, f_liq) <= ro_iwe + eps);
end
