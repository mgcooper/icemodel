function [k, dkdT] = thermal_conductivity_water(T, rho, reference)
   %THERMAL_CONDUCTIVITY_WATER Liquid-water thermal conductivity and derivative.
   %
   %  [k, dkdT] = icemodel.kernels.thermal_conductivity_water(T)
   %  [k, dkdT] = icemodel.kernels.thermal_conductivity_water(T, rho)
   %  [k, dkdT] = icemodel.kernels.thermal_conductivity_water(T, rho, reference)
   %
   % Description:
   %  Returns liquid-water thermal conductivity k [W m-1 K-1] and its
   %  constant-density temperature derivative dkdT [W m-1 K-2] for numeric
   %  arrays T [K]. The default implementation is a simple near-melting
   %  linearization of the IAPWS (2011) formulation, intended for practical
   %  snow/ice production use. The full IAPWS temperature-density
   %  formulation is also available, with the critical enhancement term
   %  omitted because this repo's application is focused on near-melting
   %  liquid water rather than near-critical water.
   %
   % Inputs:
   %  T         - Temperature [K]
   %  rho       - Liquid-water density [kg m-3]
   %  reference - Conductivity formulation (default: "iapws_2011_linear_0c")
   %
   % Options:
   %  "iapws_2011_linear_0c" - First-order Taylor expansion of the IAPWS
   %                           (2011) non-critical correlation about
   %                           T0 = 273.15 K at rho = 999.8395 kg m-3:
   %
   %                              k(T) = 0.5556444681 ...
   %                                   + 0.0024605629 * (T - 273.15)
   %
   %                           The rho input is ignored for this option.
   %  "iapws_2011" - IAPWS (2011) liquid-water conductivity using the
   %                 dilute-gas and finite-density terms only.
   %
   % Notes:
   %  The full IAPWS release defines conductivity as the sum of a
   %  temperature-density base correlation plus a critical enhancement term.
   %  For the intended snow/firn meltwater regime, the critical enhancement
   %  is negligible and is therefore omitted intentionally here.
   %
   %  The official IAPWS validity range for stable liquid water begins at the
   %  triple point. The release also states that the equation behaves in a
   %  physically reasonable manner for metastable subcooled liquid at
   %  atmospheric pressure down to at least 250 K.
   %
   % References:
   %
   %   IAPWS (2011), Release on the IAPWS Formulation 2011 for the Thermal
   %   Conductivity of Ordinary Water Substance.
   %
   % See also: thermal_conductivity_air, thermal_conductivity_ice

   arguments
      T {mustBeNumeric} = 273.15
      rho {mustBeNumeric} = 999.8395
      reference (1, 1) string {mustBeMember(reference, ...
         ["iapws_2011_linear_0c", "iapws_2011"])} = "iapws_2011_linear_0c"
   end

   switch reference
      case "iapws_2011_linear_0c"
         [k, dkdT] = iapws_2011_linear_0c(T);

      case "iapws_2011"
         [k, dkdT] = iapws_2011_noncritical(T, rho);
   end
end

function [k, dkdT] = iapws_2011_linear_0c(T)
   % First-order Taylor expansion of IAPWS 2011 about 273.15 K.

   T0 = 273.15;                       % [K]
   k0 = 0.555644468051926;            % [W m^-1 K^-1]
   dkdT0 = 0.00246056289709851;       % [W m^-1 K^-2]

   k = k0 + dkdT0 .* (T - T0);
   dkdT = dkdT0 + zeros(size(T));
end

function [k, dkdT] = iapws_2011_noncritical(T, rho)
   % IAPWS 2011 base correlation without the critical enhancement term.

   Tstar = 647.096;    % [K]
   rhostar = 322.0;    % [kg m^-3]
   kstar = 1e-3;       % [W m^-1 K^-1]

   Lk = [ ...
      2.443221e-3, ...
      1.323095e-2, ...
      6.770357e-3, ...
      -3.454586e-3, ...
      4.096266e-4];

   Lij = [ ...
      1.60397357, -0.646013523, 0.111443906, 0.102997357, -0.0504123634, ...
      0.00609859258; 2.33771842, -2.78843778, 1.53616167, -0.463045512, ...
      0.0832827019, -0.00719201245; 2.19650529, -4.54580785, 3.55777244, ...
      -1.40944978, 0.275418278, -0.0205938816; -1.21051378, 1.60812989, ...
      -0.621178141, 0.0716373224, 0.0, 0.0; -2.7203370, 4.57586331, ...
      -3.18369245, 1.1168348, -0.19268305, 0.0129138427];

   Tbar = T ./ Tstar;
   rhobar = rho ./ rhostar;

   denom = zeros(size(Tbar));
   ddenom_dTbar = zeros(size(Tbar));
   for kidx = 0:4
      denom = denom + Lk(kidx + 1) ./ Tbar .^ kidx;
      if kidx > 0
         ddenom_dTbar = ddenom_dTbar ...
            - kidx .* Lk(kidx + 1) ./ Tbar .^ (kidx + 1);
      end
   end

   lambda0_bar = sqrt(Tbar) ./ denom;
   dlambda0bar_dTbar = lambda0_bar .* (0.5 ./ Tbar - ddenom_dTbar ./ denom);

   theta = 1 ./ Tbar - 1;
   delta = rhobar - 1;

   poly = zeros(size(Tbar));
   dpoly_dTbar = zeros(size(Tbar));
   for i = 0:4
      row = zeros(size(Tbar));
      for j = 0:5
         row = row + Lij(i + 1, j + 1) .* delta .^ j;
      end

      poly = poly + theta .^ i .* row;
      if i > 0
         dpoly_dTbar = dpoly_dTbar ...
            - i .* theta .^ (i - 1) .* row ./ Tbar .^ 2;
      end
   end

   lambda1_bar = exp(rhobar .* poly);
   dlambda1bar_dTbar = lambda1_bar .* rhobar .* dpoly_dTbar;

   k = kstar .* lambda0_bar .* lambda1_bar;
   dkdT = (kstar ./ Tstar) .* ...
      (dlambda0bar_dTbar .* lambda1_bar + lambda0_bar .* dlambda1bar_dTbar);
end
