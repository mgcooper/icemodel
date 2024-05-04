function Qli = LONGIN(Tair, ea, SB)
   %LONGIN Compute incoming longwave radiation at the surface
   %
   %#codegen

   emiss_cloud = 1.08 * (1.0 - exp(-(0.01 * ea) ^ (Tair / 2016.0)));
   Qli = emiss_cloud * SB * Tair ^ 4;
end
