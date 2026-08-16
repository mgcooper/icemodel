function [radius, U_vap_faces, dm_vap] = vapor_mass_transfer(T, Ts, f_ice, f_liq, ...
      radius, dz, delz, fn, dt, varargin)
   %VAPOR_MASS_TRANSFER Compute vapor mass flux and update grain radius.
   %
   % [radius, U_vap, dm_vap] = vapor_mass_transfer(T, Ts, f_ice, f_liq, radius, ...
   %    dz, delz, fn, dt)
   % [radius, U_vap, dm_vap] = vapor_mass_transfer(..., dt, ro_vap, De)
   % [radius, U_vap, dm_vap] = ...
   %    vapor_mass_transfer(..., dt, ro_vap, De, d_vap_sfc)
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
   % Boundary conditions: zero flux at the bottom, and at the top either a
   % Dirichlet ghost node at Ts or a Neumann flux the caller supplies. See
   % the d_vap_sfc input below.
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
   %   ro_vap - (optional) Saturation vapor density at the nodes [kg m-3]
   %            (JJ x 1), from the accepted solve. Supplying it skips the
   %            exponential this function would otherwise evaluate again.
   %   De     - (optional) Effective vapor diffusivity at the nodes
   %            [m2 s-1] (JJ x 1), from the accepted solve. Supplying it
   %            skips the (T/Tf)^nd power.
   %   d_vap_sfc - (optional) Liquid-water volume fraction the surface
   %            exchanged over this step [-], positive downward, from
   %            icemodel.surface.potential_surface_vapor_exchange. Supplying
   %            it selects the Neumann top boundary described below. This
   %            function converts it to the face mass flux with
   %            icemodel.surface.surface_vapor_mass_flux, so the caller
   %            never handles the kilogram basis.
   %
   % Top boundary. Without d_vap_sfc the top face uses a Dirichlet ghost node
   % at Ts with the ice-phase saturation density. That node carries no control
   % volume, so it can supply unlimited mass, and it uses no atmospheric
   % state. With d_vap_sfc the top face carries the turbulent surface exchange
   % the surface energy balance already computed, and no ghost node is
   % evaluated. That is the boundary the coupled vapor mode uses: the surface
   % exchange is computed once, in the SEB, and applied once, here.
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

   persistent nd g1 g2 r_max Uv_max
   if isempty(nd)
      [nd, g1, g2, r_max, Uv_max] = icemodel.parameterLookup( ...
         'nd', 'g1', 'g2', 'r_max', 'Uv_max');
   end

   JJ = numel(T);

   % The caller supplies the accepted solve's node quantities, or this
   % function evaluates them. Reuse is what keeps the exponential and the
   % power off the substep path a second time.
   use_neumann_top = nargin > 11;

   % --- Saturation vapor density at each node ---

   % Phase-aware vapor density [kg m-3] and diffusivity [m2 s-1] at each node.
   if nargin > 9
      ro_vap = varargin{1};
   else
      ro_vap = icemodel.vapor.saturation_vapor_density(T, f_liq);
   end
   if nargin > 10
      De = varargin{2};
   else
      De = icemodel.vapor.vapor_diffusivity(T);
   end

   % --- Top boundary ---

   if use_neumann_top
      % The surface exchange arrives as a face flux, so pad with node 1 and
      % overwrite face 1 below. This evaluates no ghost node, which is the
      % point: the saturated ghost can supply unlimited mass and reads no
      % atmospheric state.
      ro_vap_s = ro_vap(1);
      De_s = De(1);
   else
      % Dirichlet ghost node. Use ice-phase es at surface (sublimating
      % interface).
      f_liq_s = 0;
      ro_vap_s = icemodel.vapor.saturation_vapor_density(Ts, f_liq_s);
      De_s = icemodel.vapor.vapor_diffusivity(Ts);
   end

   % --- Vapor flux at control volume interfaces (Patankar Eq. 4.9) ---

   % Vapor mass flux at the JJ+1 interfaces [kg m-2 s-1].
   % Positive = downward (from surface into column). The bottom face is zero
   % flux (Neumann / insulated). Face diffusivity is the fn-weighted harmonic
   % mean of the node values, Patankar Eq. 4.9.
   U_vap_faces = icemodel.column.vapor_face_quantities( ...
      ro_vap, ro_vap_s, De, De_s, delz, fn);

   % The turbulent exchange replaces the diffusive top face. It is applied
   % once, here, having been computed once, in the surface energy balance.
   % The caller passes a fraction, so the kilogram basis is formed here.
   if use_neumann_top
      U_vap_faces(1) = icemodel.surface.surface_vapor_mass_flux( ...
         varargin{3}, dz(1), dt);
   end

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
