function Qa = advective_heat_flux(ppt, tppt, cv_liq)
   %ADVECTIVE_HEAT_FLUX Compute heat advected to the surface by rainfall.
   %
   %  Qa = icemodel.surface.advective_heat_flux(ppt, tppt, cv_liq)
   %
   % Computes the heat flux [W m-2] due to rainfall at rate ppt [m s^-1]
   % and wet-bulb temperature tppt [K]:
   %
   %   Qa = cv_liq * ppt * tppt
   %
   % where cv_liq = cp_liq * ro_liq [J m^-3 K^-1].
   %
   % This formula is valid for rain only (Jordan 1991, Eq. A.4). Snowfall
   % requires a separate treatment using cp_snow * ro_snow * fallrate * tppt.
   % ppt must be in m s^-1 (a rainfall rate, not an accumulation depth).
   %
   %#codegen
   Qa = cv_liq * ppt * tppt;
end
