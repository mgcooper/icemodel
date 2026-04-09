function Qa = advective_heat_flux(ppt, tppt, cv_liq)
   %ADVECTIVE_HEAT_FLUX Compute heat due to advection of liquid water.
   %
   %  Qa = icemodel.surface.advective_heat_flux(ppt, tppt, cv_liq)
   %
   % Computes heat flux [W m-2] due to liquid water with precipitation rate
   % ppt [m] at temperature tppt [K] advected onto the surface.
   %
   %#codegen
   Qa = cv_liq * ppt * tppt;
end
