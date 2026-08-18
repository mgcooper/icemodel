function [radius, U_vap_faces, dm_vap] = vapor_mass_transfer(T, Ts, f_ice, f_liq, ...
      radius, dz, delz, fn, dt, varargin)
   %VAPOR_MASS_TRANSFER Compute vapor mass flux and update grain radius.
   %
   % [radius, U_vap, dm_vap] = vapor_mass_transfer(T, Ts, f_ice, f_liq, radius, ...
   %    dz, delz, fn, dt)
   % [radius, U_vap, dm_vap] = vapor_mass_transfer(..., dt, d_vap_faces)
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
   % Boundary conditions: zero flux at the bottom and a Dirichlet ghost
   % node at Ts on top. The coupled path below computes no fluxes at all:
   % it consumes the accumulated face exchange the driver supplies.
   %
   % Grain growth follows Jordan (1991) in terms of diameter d = 2*r:
   %   Dry snow (f_liq < 1e-4):   dd/dt = g1 / d * |U_vap|    (Eq. 33)
   %   Wet, low (f_liq < 0.09):   dd/dt = g2 / d * (f_liq + 0.05) (Eq. 34a)
   %   Wet, high (f_liq >= 0.09): dd/dt = g2 / d * 0.14        (Eq. 34b)
   %
   % The Jordan grain growth model is monotonic. Grains grow and never shrink.
   % Shrinkage needs other mechanisms that this function does not include:
   % fresh snow deposition that resets the grain size, wind slab formation, or
   % surface rounding driven by sublimation.
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
   %   d_vap_faces - (optional) Gross face exchange the accepted substeps
   %            applied over this step [m w.e.] (JJ+1 x 1): the realized
   %            surface exchange at face 1 and the substep-integrated
   %            interior transport magnitudes elsewhere, from the driver's
   %            accumulation through icemodel.column.couple_vapor_step.
   %            Supplying it selects the coupled grain-growth path below.
   %
   % Coupled path. With d_vap_faces this function grows grains from the
   % fluxes the column actually transported, per DesignSpec decision 8: the
   % accumulated depths convert to step-mean magnitude fluxes, and no
   % saturation state, ghost node, or face flux is evaluated at all, so the
   % path adds no exponential and no power. The magnitudes are gross so
   % substeps whose exchange reversed sign add rather than cancel, and the
   % surface entry is the realized exchange, so demand the applier rejected
   % never drives growth. U_vap in this mode returns those step-mean
   % magnitudes, and dm_vap returns zeros: magnitudes carry no sign, so no
   % divergence exists on this path.
   %
   % Default path. Without d_vap_faces the fluxes are computed here: the top
   % face uses a Dirichlet ghost node at Ts with the ice-phase saturation
   % density. That node carries no control volume, so it can supply
   % unlimited mass, and it uses no atmospheric state.
   %
   % Outputs:
   %   radius - Updated grain effective radius [m] (JJ x 1)
   %   U_vap  - Vapor mass flux at interfaces [kg m-2 s-1] (JJ+1 x 1)
   %              Positive = downward (into column from surface).
   %   dm_vap - Volumetric mass source rate [kg m-3 s-1] (JJ x 1)
   %              Positive = deposition (mass gain), negative = sublimation.
   %
   % This function uses radius to match the spectral model, which uses the
   % optically equivalent grain radius. The thermal grain radius tracked here
   % and the spectral radius are not the same quantity. Coupling them is
   % future work. See icemodel.radiation.initialize_spectral_model,
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

   persistent nd g1 g2 r_max Uv_max ro_liq
   if isempty(nd)
      [nd, g1, g2, r_max, Uv_max] = icemodel.parameterLookup( ...
         'nd', 'g1', 'g2', 'r_max', 'Uv_max');
      ro_liq = icemodel.physicalConstant('ro_liq');
   end

   JJ = numel(T);

   % The driver supplies the accumulated step exchange, or this function
   % computes its own fluxes. The accumulation is what lets the coupled
   % path skip every saturation evaluation.
   use_step_fluxes = nargin > 9;

   if use_step_fluxes

      % --- Coupled path: consume the accumulated face exchange ---

      % Convert the gross water-equivalent depths to step-mean magnitude
      % fluxes [kg m-2 s-1]. These are the fluxes the substeps applied:
      % realized at the surface, transported in the interior.
      U_vap_faces = varargin{1} * ro_liq / dt;

      % Magnitudes carry no sign, so no flux divergence exists on this path.
      dm_vap = zeros(JJ, 1);
   else

      % --- Saturation vapor density at each node ---

      % Phase-aware vapor density [kg m-3] and diffusivity [m2 s-1] at each
      % node.
      ro_vap = icemodel.vapor.saturation_vapor_density(T, f_liq);
      De = icemodel.vapor.vapor_diffusivity(T);

      % --- Top boundary ---

      % Dirichlet ghost node. Use ice-phase es at surface (sublimating
      % interface).
      f_liq_s = 0;
      ro_vap_s = icemodel.vapor.saturation_vapor_density(Ts, f_liq_s);
      De_s = icemodel.vapor.vapor_diffusivity(Ts);

      % --- Vapor flux at control volume interfaces (Patankar Eq. 4.9) ---

      % Vapor mass flux at the JJ+1 interfaces [kg m-2 s-1].
      % Positive = downward (from surface into column). The bottom face is
      % zero flux (Neumann / insulated). Face diffusivity is the fn-weighted
      % harmonic mean of the node values, Patankar Eq. 4.9.
      U_vap_faces = icemodel.column.vapor_face_quantities( ...
         ro_vap, ro_vap_s, De, De_s, delz, fn);

      % --- Mass source from flux divergence ---

      % Net flux into each control volume [kg m-3 s-1]
      dm_vap = (U_vap_faces(1:JJ) - U_vap_faces(2:JJ+1)) ./ dz;
   end

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
