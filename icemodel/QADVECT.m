function Qa = QADVECT(ppt,tppt,cv_liq)
   %QADVECT Compute heat due to advection of liquid water
   %
   %  Qa = QADVECT(ppt, tppt, cv_liq) Computes heat flux [W m-2] due to liquid
   %  water (ppt) [m] with temperature [K] advected onto the surface.
   Qa = cv_liq * ppt * tppt;
end
