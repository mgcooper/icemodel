function nu_air = air_kinematic_viscosity(Ta, ro_atm)
   %AIR_KINEMATIC_VISCOSITY Approximate air kinematic viscosity.
   %
   % Dynamic viscosity is approximated with a Sutherland-law form:
   %   mu_air = mu0 * (Ta / T0)^(3/2) * (T0 + S) / (Ta + S)
   %
   % where:
   %   mu0    = reference dynamic viscosity [Pa s]
   %   T0     = reference temperature [K]
   %   S      = Sutherland constant [K]
   %   Ta     = air temperature [K]
   %
   % Kinematic viscosity then follows from:
   %   nu_air = mu_air / ro_atm
   %
   % where ro_atm is the local moist-air density [kg m^-3].
   %
   %#codegen

   mu0 = 1.716e-5;
   T0 = 273.15;
   S = 111.0;

   mu_air = mu0 * (Ta / T0) ^ (3 / 2) * (T0 + S) / (Ta + S);
   nu_air = mu_air / ro_atm;
end
