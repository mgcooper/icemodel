function [q, dq_df_liq] = liquid_flux(f_liq, f_ice, kwargs)
   %LIQUID_FLUX Compute the liquid water flux between snowpack layers.
   %
   %  q                = icemodel.column.liquid_flux(f_liq, f_ice)
   %  [q, dq_df_liq]   = icemodel.column.liquid_flux(f_liq, f_ice, ...)
   %
   %  Returns the volumetric liquid water flux q [m s-1] (per unit snow column
   %  area, integrated over a substep at the call site by multiplying by dt) and
   %  its derivative dq/d(f_liq) [s-1] for the supplied column state. The
   %  constitutive law evaluated is the Mualem-style power law on the
   %  colbeck1972 / darcy paths
   %
   %     q = k_sat * relSat ^ m_exp
   %
   %  and the van-Genuchten retention form on the shimizu1970 path. The exponent
   %  m_exp is loaded from icemodel.parameterLookup ("m_exp", default 3) per
   %  Clark et al. (2017) Eq. 7 and Clark et al. (2021) Table 4.
   %
   %  Parameters:
   %  -----------
   %  f_liq : vector
   %      Volumetric liquid water fraction in each layer [1].
   %
   %  f_ice : vector
   %      Volumetric ice fraction in each layer [1].
   %
   %  Name-value:
   %  -----------
   %  f_res_pore : scalar (optional, also referred to as liqresid or Swi)
   %      Residual water volume per pore volume [1]. Production callers pass
   %      f_res_pore_snow/ice/firn via icemodel.surface.update_surface_state.
   %      Default NaN routes to icemodel.parameterLookup('f_res_pore'), the
   %      Colbeck-benchmark value (0.07).
   %
   %  k_sat_method : string (optional, default "colbeck1972")
   %      "colbeck1972" | "shimizu1970" | "darcy". Selects the saturated
   %      hydraulic conductivity closure scheme; the q-vs-S retention form is
   %      power-law for "colbeck1972" / "darcy" and van-Genuchten for
   %      "shimizu1970".
   %
   %  grainsz : scalar (optional)
   %      Grain diameter [m]; required when k_sat_method is "shimizu1970".
   %
   %  permeability : scalar (optional)
   %      Intrinsic permeability kappa [m^2]; required when k_sat_method is
   %      "darcy". Sets k_sat = kappa * ro_liq * g / mu_water_ref.
   %
   %  Returns:
   %  --------
   %  q : vector
   %      Liquid water flux [m s-1] between the layers.
   %
   %  dq_df_liq : vector
   %      Derivative dq/d(f_liq) [s-1].
   %
   %  Notes:
   %  ------
   %  This function dispatches between three saturated-hydraulic- conductivity
   %  closures:
   %    1. Colbeck 1972 (density-only, default)
   %    2. Shimizu 1970 (grain-size + water-fraction; van-Genuchten q-S)
   %    3. Darcy form  k_sat = kappa * ro_liq * g / mu  (caller-pinned
   %                   permeability; matches the SUMMA Laugh-Tests
   %                   parameter trial values for the Colbeck case)
   %
   %  The function calculates available capacity and relative saturation of
   %  liquid water in the snowpack. The flux is then computed from the saturated
   %  hydraulic conductivity and relative saturation.
   %
   %  Variable analogues with Colbeck 1972 (used throughout this file):
   %
   %      f_res_pore - residual water volume / pore volume [1] (Swi)
   %      f_por      - pore fraction [1]                       (phi)
   %      f_res      - residual water fraction [1]
   %      liqsat     - water saturation [1]                    (Sw)
   %      availCap   - available capacity [1]
   %      relSat     - relative saturation [1]                 (Sstar)
   %
   %  f_ice, f_liq, and f_air are volumetric fractions. f_wat is the
   %  water-equivalent volume fraction, i.e., if the ice melted, the combined
   %  liquid plus melted-ice fraction.
   %
   %  Note: in general, 2-7 % of the pore space must be filled with water before
   %  any can infiltrate, so when debugging, check that the residual capillary
   %  floor is not exceeded.
   %
   %  References
   %    Colbeck, S. C. (1972). A theory of water percolation in snow.
   %    Shimizu, H. (1970). Air permeability of deposited snow.
   %    Clark, M. P. et al. (2017). An analytical test case for snow
   %      models.
   %    Clark, M. P. et al. (2021). The numerical implementation of
   %      land models: Problem formulation and laugh tests.
   %
   % See also: icemodel.column.saturated_hydraulic_conductivity,
   %  icemodel.column.infiltration, icemodel.parameterLookup
   %
   %#codegen

   arguments
      f_liq
      f_ice
      kwargs.f_res_pore     (1, 1) double = NaN
      kwargs.k_sat_method   (1, 1) string {mustBeMember(kwargs.k_sat_method, ...
         ["colbeck1972", "shimizu1970", "darcy"])} = "colbeck1972"
      kwargs.grainsz        (1, 1) double = NaN
      kwargs.permeability   (1, 1) double = NaN
   end

   % Load constitutive-law constants once per session. m_exp is the Mualem
   % exponent (Clark 2017 Eq. 7, Clark 2021 Table 4); shi_n_a / shi_n_b are the
   % van-Genuchten n shape parameters (slope and grainsz exponent); shi_a_wat is
   % the Shimizu water-fraction exponent that also enters the closed-form
   % dk_sat/df_liq term.
   persistent m_exp default_f_res_pore shi_n_a shi_n_b shi_a_wat
   if isempty(m_exp)
      [m_exp, default_f_res_pore, shi_n_a, shi_n_b, shi_a_wat] = ...
         icemodel.parameterLookup('m_exp', 'f_res_pore', 'shi_n_a', ...
         'shi_n_b', 'shi_a_wat');
   end

   % Resolve the residual capillary-saturation kwarg. NaN is the
   % verification-default sentinel: production callers pass an explicit
   % per-substep value derived from opts.f_res_pore_snow / _ice / _firn.
   if isnan(kwargs.f_res_pore)
      f_res_pore = default_f_res_pore;
   else
      f_res_pore = kwargs.f_res_pore;
   end

   % Initialize outputs to the no-flow state. Layer indices that satisfy the
   % flow criterion below are then filled in place; layers at or below the
   % residual capillary floor remain at q = 0.
   N = numel(f_liq);
   q = zeros(N, 1);
   dq_df_liq = zeros(N, 1);

   % Saturated hydraulic conductivity from the column-level kernel.
   % This dispatches on k_sat_method to one of:
   %   colbeck1972 - density-only fit  k_sat = col_k0 * exp(col_a_phi * f_por)
   %   shimizu1970 - grain-size + f_wat fit (varies in time as f_wat moves)
   %   darcy       - constant k_sat from caller-supplied permeability
   %
   % If updating dynamic viscosity as a function of T:
   %   n = icemodel.dynamicViscosityWater(T);
   %   k_sat = (6.1313e-10 ./ n) .* exp(15.9 * f_por);
   k_sat = icemodel.column.saturated_hydraulic_conductivity( ...
      f_ice, f_liq, ...
      method=kwargs.k_sat_method, ...
      grainsz=kwargs.grainsz, ...
      permeability=kwargs.permeability);

   % Pore fraction, residual liquid fraction, available pore capacity, and
   % relative saturation. These are the same SUMMA-style quantities used
   % elsewhere in icemodel.column. relSat is set to 1 in fully-ice layers
   % (availCap == 0) so the no-flux mask below picks them up cleanly.
   f_por = 1.0 - f_ice;
   f_res = f_res_pore * f_por;
   availCap = max(0.0, f_por - f_res);
   relSat = (f_liq - f_res) ./ availCap;
   relSat(availCap <= 0) = 1;

   % Note: if relSat is used elsewhere, set negative values to 0 or NaN:
   % relSat(relSat < 0) = NaN;

   % Indices where flow can occur. Layers below the residual capillary floor or
   % with no pore space contribute zero flux.
   iflux = availCap > 0 & f_liq > f_res;

   switch kwargs.k_sat_method

      case {"colbeck1972", "darcy"}
         % Power-law q-vs-S retention. Both methods share this branch because
         % they differ only in how k_sat is computed (formula vs caller-pinned
         % permeability); the constitutive q(S) and its derivative are
         % identical.
         q(iflux) = k_sat(iflux) .* relSat(iflux) .^ m_exp;

         % Closed-form derivative dq/d(f_liq):
         %   dq/dS = m_exp * k_sat * S^(m_exp - 1)
         %   dS/d(f_liq) = 1 / availCap
         dq_df_liq(iflux) = m_exp .* k_sat(iflux) ...
            .* relSat(iflux) .^ (m_exp - 1) ./ availCap(iflux);

      case "shimizu1970"
         % Van-Genuchten q-vs-S form coupled to the Shimizu permeability. Both
         % share grainsz and are kept tied as one method here.
         if isnan(kwargs.grainsz)
            error('icemodel:column:liquid_flux:missingGrainsz', ...
               'shimizu1970 method requires the grainsz keyword input [m].');
         end

         % Van-Genuchten n and m shape parameters. n is grain-size dependent
         % (shi_n_a * exp(-shi_n_b * 2 * grainsz) + 1).
         n_vg = shi_n_a * exp(-shi_n_b * 2 * kwargs.grainsz) + 1;
         m_vg = 1 - 1 / n_vg;

         % q(S) = k_sat * sqrt(S) * (1 - (1 - S^(1/m_vg))^m_vg)^2
         S = relSat(iflux);
         A = 1 - S .^ (1 / m_vg);
         B = 1 - A .^ m_vg;
         q(iflux) = k_sat(iflux) .* sqrt(S) .* B .^ 2;

         % Closed-form derivative dq/d(f_liq). On the Shimizu path k_sat depends
         % on f_liq via f_wat:
         %   d k_sat / d f_liq = -shi_a_wat * k_sat
         % and the saturation chain rule gives
         %   dq/dS    = 0.5 / sqrt(S) * B^2 + sqrt(S) * d(B^2)/dS
         %   d(B^2)/dS = 2 B dB/dS = 2 B * A^(m_vg-1) * S^(1/m_vg - 1)
         %   dS / d f_liq = 1 / availCap
         dB2_dS = 2 .* B .* A .^ (m_vg - 1) .* S .^ (1 / m_vg - 1);
         dq_dS = k_sat(iflux) .* (0.5 ./ sqrt(S) .* B .^ 2 + sqrt(S) .* dB2_dS);
         dq_df_liq(iflux) = -shi_a_wat .* q(iflux) + dq_dS ./ availCap(iflux);
   end
end

% For reference, in terms of CV thicknesses:
% Swi = f_res_por = f_res/f_por  = h_res/h_por = h_res/(h_tot - h_ice)
% Sw  = liqsat    = f_liq/f_por  = h_liq/h_por = h_liq/(h_tot - h_ice)
% phi = porosity  = f_por        = h_por/h_tot = (h_tot - h_ice)/h_tot
%
% Sstar = relSat  = (f_liq - f_res)/availCap = (h_liq - h_res)/availCap
