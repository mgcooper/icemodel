function [radius, U_vap_faces, dm_vap] = VAPORTRANSFER(T, Ts, f_ice, f_liq, ...
      radius, dz, delz, fn, dt)
   %VAPORTRANSFER Compute vapor mass flux and update grain radius.
   %
   % [radius, U_vap, dm_vap] = VAPORTRANSFER(T, Ts, f_ice, f_liq, radius, ...
   %    dz, delz, fn, dt)
   %
   % Computes the diffusive water vapor mass flux through the porous ice
   % column, the mass source/sink from flux divergence, and updates grain
   % radius following Jordan (1991) SNTHERM89 Eqs. 20-21, 33-34.
   %
   % The vapor flux follows Fick's law applied to the saturation vapor density:
   %
   %   U_vap = -De * d(ro_vap)/dz   [kg m-2 s-1]
   %
   % where De = De0 * (T / Tf)^nd is the effective vapor diffusivity (Yen's
   % enhancement) and ro_vap = es / (Rv * T) is the saturation vapor density
   % from the ideal gas law. Saturation vapor pressure es is obtained from
   % icemodel.vapor.saturation_vapor_pressure (Ambaum 2020 / Romps 2021
   % Rankine-Kirchhoff formula).
   %
   % Interface diffusivities use Patankar (1980) harmonic mean (Eq. 4.9).
   % Boundary conditions: Dirichlet (surface Ts) at top, zero-flux at bottom.
   %
   % Grain growth follows Jordan (1991) in terms of diameter d = 2*r:
   %   Dry snow (f_liq < 1e-4):   dd/dt = g1 / d * |U_vap|    (Eq. 33)
   %   Wet, low (f_liq < 0.09):   dd/dt = g2 / d * (f_liq + 0.05) (Eq. 34a)
   %   Wet, high (f_liq >= 0.09): dd/dt = g2 / d * 0.14        (Eq. 34b)
   %
   % Note: the Jordan grain growth model is strictly monotonic — grains grow
   % but never shrink. Physically, shrinkage requires additional mechanisms
   % (e.g., fresh snow deposition resetting grain size, wind slab formation,
   % or sublimation-driven surface rounding) that are not included here.
   %
   % Inputs:
   %   T      - Node temperatures [K] (JJ x 1)
   %   Ts     - Surface temperature [K] (scalar)
   %   f_ice  - Volumetric ice fraction (JJ x 1)
   %   f_liq  - Volumetric liquid water fraction (JJ x 1)
   %   radius - Grain effective radius [m] (JJ x 1)
   %   dz     - Control volume thicknesses [m] (JJ x 1)
   %   delz   - Distances between adjacent node centers [m] (JJ+1 x 1)
   %   fn     - Interface interpolation weights (JJ+1 x 1)
   %   dt     - Timestep [s]
   %
   % Outputs:
   %   radius - Updated grain effective radius [m] (JJ x 1)
   %   U_vap  - Vapor mass flux at interfaces [kg m-2 s-1] (JJ+1 x 1)
   %              Positive = downward (into column from surface).
   %   dm_vap - Volumetric mass source rate [kg m-3 s-1] (JJ x 1)
   %              Positive = deposition (mass gain), negative = sublimation.
   %
   % Note: radius is used for consistency with the spectral model (optically
   % equivalent grain radius). Long-term, the thermal grain radius tracked
   % here and the spectral radius are not identical quantities; coupling them
   % is future work. See icemodel.radiation.initialize_spectral_model,
   % icemodel.radiation.update_extinction_coefficients.
   %
   % References:
   %   Jordan (1991), "A one-dimensional temperature model for a snow cover:
   %      Technical documentation for SNTHERM.89." CRREL Special Report 91-16.
   %   Patankar (1980), "Numerical Heat Transfer and Fluid Flow." CRC Press.
   %
   % See also: icemodel.vapor.saturation_vapor_density,
   %  icemodel.vapor.vapor_diffusivity,
   %  icemodel.vapor.saturation_vapor_pressure,
   %  icemodel.column.bulk_thermal_conductivity,
   %  icemodel.radiation.initialize_spectral_model
   %
   %#codegen

   persistent nd g1 g2 r_max Uv_max
   if isempty(nd)
      [nd, g1, g2, r_max, Uv_max] = icemodel.parameterLookup( ...
         'nd', 'g1', 'g2', 'r_max', 'Uv_max');
   end

   JJ = numel(T);

   % --- Saturation vapor density at each node ---

   % Phase-aware vapor density [kg m-3] and diffusivity [m2 s-1] at each node.
   ro_vap = icemodel.vapor.saturation_vapor_density(T, f_liq);
   De = icemodel.vapor.vapor_diffusivity(T);

   % --- Surface ghost node ---

   % Use ice-phase es at surface (sublimating interface)
   f_liq_s = 0;
   ro_vap_s = icemodel.vapor.saturation_vapor_density(Ts, f_liq_s);
   De_s = icemodel.vapor.vapor_diffusivity(Ts);

   % --- Vapor flux at control volume interfaces (Patankar Eq. 4.9) ---

   % Padded arrays: [surface; nodes 1:JJ; bottom ghost]
   ro_vap_nodes = [ro_vap_s; ro_vap; ro_vap(JJ)];
   De_nodes = [De_s; De; De(JJ)];

   % Harmonic mean diffusivity at JJ+1 interfaces
   De_faces = 1.0 ./ ...
      ((1.0 - fn) ./ De_nodes(1:JJ+1) + fn ./ De_nodes(2:JJ+2));

   % Vapor mass flux at interfaces [kg m-2 s-1]
   % Positive = downward (from surface into column)
   U_vap_faces = -De_faces .* ...
      (ro_vap_nodes(2:JJ+2) - ro_vap_nodes(1:JJ+1)) ./ delz;

   % Bottom boundary: zero flux (Neumann / insulated)
   U_vap_faces(JJ+1) = 0;

   % --- Mass source from flux divergence ---

   % Net flux into each control volume [kg m-3 s-1]
   dm_vap = (U_vap_faces(1:JJ) - U_vap_faces(2:JJ+1)) ./ dz;

   % --- Grain growth (Jordan 1991, SNTHERM89 Eqs. 33-34) ---

   % Jordan equations are in diameter d = 2*r. Convert at boundaries.
   diam = 2 * radius;

   % Vapor flux magnitude at nodes: average of adjacent interface magnitudes
   U_vap_nodes = min( ...
      0.5 * (abs(U_vap_faces(1:JJ)) + abs(U_vap_faces(2:JJ+1))), Uv_max);

   % Dry snow: vapor-driven growth (Eq. 33)
   dry = f_liq < 1e-4;
   diam(dry) = diam(dry) + dt * g1 .* U_vap_nodes(dry) ./ diam(dry);

   % Wet snow, moderate liquid (Eq. 34a)
   wet_lo = ~dry & f_liq < 0.09;
   diam(wet_lo) = ...
      diam(wet_lo) + dt * g2 .* (f_liq(wet_lo) + 0.05) ./ diam(wet_lo);

   % Wet snow, high liquid (Eq. 34b)
   wet_hi = f_liq >= 0.09;
   diam(wet_hi) = diam(wet_hi) + dt * g2 * 0.14 ./ diam(wet_hi);

   % Convert back to radius, enforce maximum
   radius = min(0.5 * diam, r_max);
end
