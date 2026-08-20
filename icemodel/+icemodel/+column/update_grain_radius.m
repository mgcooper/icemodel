function r_eff = update_grain_radius(r_eff, f_liq, U_vap, d_vap_applied, ...
      dz1, dt)
   %UPDATE_GRAIN_RADIUS Grow thermal grains from the substep vapor exchange.
   %
   % r_eff = update_grain_radius(r_eff, f_liq, U_vap, d_vap_applied, dz1, dt)
   %
   % Updates the thermal grain radius following Jordan (1991) SNTHERM89
   % Eqs. 33-34, once per accepted substep. Dry growth uses the substep's
   % own instantaneous face-flux magnitudes, the form Jordan's Eq. 33
   % integrates, so no cross-substep accumulator is needed. Wet growth
   % uses the concurrent liquid fraction. This function does not diagnose
   % vapor transport or saturation state.
   %
   % Grain growth is monotonic. Other processes must represent grain-size
   % reductions, including fresh-snow deposition, wind-slab formation, and
   % sublimation-driven surface rounding.
   %
   % Inputs:
   %   r_eff         - Thermal grain effective radius [m] (JJ x 1)
   %   f_liq         - Volumetric liquid water fraction (JJ x 1)
   %   U_vap         - Accepted interior vapor mass flux [kg m-2 s-1]
   %                   (JJ+1 x 1, boundary faces zero)
   %   d_vap_applied - Realized surface vapor exchange this substep [-],
   %                   a signed top-cell liquid-water volume fraction
   %   dz1           - Top control-volume thickness [m]
   %   dt            - Accepted substep duration [s]
   %
   % Output:
   %   r_eff         - Updated thermal grain effective radius [m] (JJ x 1)
   %
   % The thermal radius is distinct from the optically equivalent radius in
   % the spectral model. See icemodel.radiation.initialize_spectral_model and
   % icemodel.radiation.update_extinction_coefficients.
   %
   % Reference:
   %   Jordan (1991), "A one-dimensional temperature model for a snow cover:
   %      Technical documentation for SNTHERM.89." CRREL Special Report 91-16.
   %
   % See also: icemodel.column.couple_vapor_step
   %
   %#codegen

   persistent g1 g2 r_max Uv_max ro_liq
   if isempty(g1)
      [g1, g2, r_max, Uv_max] = icemodel.parameterLookup( ...
         'g1', 'g2', 'r_max', 'Uv_max');
      ro_liq = icemodel.physicalConstant('ro_liq');
   end

   % Growth scales with the flux magnitude, so the surface face carries the
   % realized surface exchange as a magnitude flux: rejected demand never
   % crossed the surface, and a sign never cancels interior transport.
   JJ = numel(r_eff);
   U_vap_faces = U_vap;
   U_vap_faces(1) = abs(d_vap_applied) * dz1 * ro_liq / dt;
   U_vap_nodes = min( ...
      0.5 * (abs(U_vap_faces(1:JJ)) + abs(U_vap_faces(2:JJ+1))), Uv_max);

   % Jordan's equations use grain diameter rather than radius.
   diam = 2 * r_eff;

   % Dry snow grows from realized vapor throughput (Jordan Eq. 33).
   dry = f_liq < 1e-4;
   diam(dry) = diam(dry) + dt * g1 .* U_vap_nodes(dry) ./ diam(dry);

   % Moderately wet snow grows at a liquid-dependent rate (Jordan Eq. 34a).
   wet_lo = ~dry & f_liq < 0.09;
   diam(wet_lo) = ...
      diam(wet_lo) + dt * g2 .* (f_liq(wet_lo) + 0.05) ./ diam(wet_lo);

   % Wetter snow uses the capped liquid-growth rate (Jordan Eq. 34b).
   wet_hi = f_liq >= 0.09;
   diam(wet_hi) = diam(wet_hi) + dt * g2 * 0.14 ./ diam(wet_hi);

   % Return radius and enforce Jordan's maximum thermal grain size.
   r_eff = min(0.5 * diam, r_max);
end
