function k = thermal_conductivity_ice(T, reference)
   %THERMAL_CONDUCTIVITY_ICE Compute thermal conductivity of ice from temperature.
   %
   %  K = icemodel.kernels.thermal_conductivity_ice(T)
   %  K = icemodel.kernels.thermal_conductivity_ice(T, reference)
   %
   % Description:
   %  Returns the thermal conductivity K [W m-1 K-1] of ice as a function of
   %  temperature T [K]. The supported options include near-melting ice-Ih
   %  fits, simple inverse-temperature approximations for pure ice Ih, and
   %  the Yen (1981) glacier/polycrystalline fit used elsewhere in this repo.
   %
   % Inputs:
   %  T         - Temperature [K] (default: 273.15)
   %  reference - Empirical reference to use (default: "andersson_2005")
   %
   % Options:
   %  "andersson_2005" - Andersson and Inaba (2005) linear fit for ice Ih
   %                     near melting.
   %  "yen_1981"   - Yen (1981), fit for glacier/polycrystalline ice.
   %                 Valid range: ~200-270 K.
   %  "rabin_2000" - Rabin (2000), power-law fit for pure water ice (Ih).
   %                 Valid range: ~100-270 K.
   %  "petrenko_1999" - Petrenko and Whitworth (1999) inverse-temperature
   %                    approximation for pure water ice (Ih).
   %  "slack_1980" - Slack (1980) / Petrenko & Whitworth (1999), the standard
   %                 inverse-temperature approximation for pure water ice
   %                 (Ih). Valid range: ~150-270 K.
   %  "engineering_toolbox" - Approximate handbook-style fit retained only
   %                          for comparison with tabulated engineering values.
   %
   % References:
   %
   %   Yen, Y.C. (1981), Review of thermal properties of snow, ice, and sea ice,
   %   CRREL Report 81-10, U.S. Army Cold Regions Research and Engineering Lab.
   %
   %   Ritz, C. (1987), Time dependent boundary conditions for calculation of
   %   temperature fields in ice sheets, in The Physical Basis of Ice Sheet
   %   Modelling, vol. 170, pp. 207-216, Vancouver.
   %
   %  Pure water ice (Ih):
   %   Rabin, Y. (2000), The effect of temperature-dependent thermal
   %   conductivity in heat transfer simulations of frozen biomaterials.
   %   https://www.andrew.cmu.edu/user/yr25/RabinPublications/Rabin_Pub230.pdf
   %
   %   Slack, G.A. (1980), Thermal conductivity of ice, Phys. Rev. B, 22(6),
   %   3065-3071.
   %
   %   Harvey (2018), Properties of Ice and Supercooled Water (NIST).
   %
   %   Petrenko, V.F. and Whitworth, R.W. (1999), Physics of Ice, Oxford
   %   University Press. (Standard reference; tabulates k vs T for ice Ih.)
   %
   % See also: thermal_conductivity_air, thermal_conductivity_snow,
   % thermal_conductivity_water

   arguments (Input)
      T {mustBeNumeric} = 273.15
      reference (1, 1) string {mustBeMember(reference, ...
         ["andersson_2005", "yen_1981", "rabin_2000", "petrenko_1999", ...
         "slack_1980", "engineering_toolbox"])} = "andersson_2005"
   end

   switch reference

      case "andersson_2005"
         % Andersson & Inaba (2005) linear fit for pure ice Ih near melting.
         % Used by Harvey (2018, NIST) table values.
         % k(273.15 K) = 2.17 W m^-1 K^-1
         %
         % Reference:
         % Andersson, O. and Inaba, A. (2005),
         % Thermal conductivity of crystalline and amorphous ices.
         k = 2.1725 - 0.003403 * (T - 273.15);

      case "yen_1981"
         % Yen (1981) exponential fit for polycrystalline ice.
         % k(273.15) = 2.07

         k = 9.828 * exp(-0.0057 * T);

      case "rabin_2000"
         % Rabin et al. (2000) power-law fit for pure ice Ih based on
         % laboratory conductivity data from ~80-220 K.
         % k(273.15) = 2.09 W m^-1 K^-1

         k = 2135 * T .^ -1.235;

      case "petrenko_1999"
         % Petrenko & Whitworth (1999) inverse-temperature approximation
         % for pure ice Ih.
         % Simple phonon-limited conductivity form:
         % k = A / T with A = 651 W m^-1
         % k(273.15) = 2.38

         k = 651 ./ T;

      case "slack_1980"
         % Slack (1980) / NIST-style fit for pure ice Ih at ~1 atm.
         % Recommended for pure frozen water near environmental temperatures.
         % k(273.15) = 2.26 W m^-1 K^-1

         k = 488.19 ./ T + 0.4685;

      case "engineering_toolbox"
         % Power-law fit through standard engineering handbook tabulated
         % ice-Ih values, anchored at k(273.15) = 2.24 W m^-1 K^-1.
         % Useful only if consistency with handbook tables is desired.
         %
         % Approximate source:
         % Engineering Toolbox tabulated ice thermal conductivity values.
         k = 2.24 * (273.15 ./ T) .^ 1.885;
   end
end
