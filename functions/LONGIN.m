function Qli = LONGIN(ea,Tair,Stef_Boltz)
   %LONGIN Compute incoming longwave radiation at the surface

   emiss_cloud = 1.08 * (1.0 - exp(-(0.01 * ea)^(Tair/2016.)));
   Qli = emiss_cloud * Stef_Boltz * Tair^4;
end
