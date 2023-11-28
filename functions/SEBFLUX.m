function [Qe, Qh, Qc, Qm, Qf, balance] = SEBFLUX(T, xTsfc, Tair, Qsi, Qli, ...
      albedo, wspd, Pa, De, ea, Tf, k_eff, dz, cv_air, roL, emiss, SB, chi, ...
      epsilon, scoef, liqflag) %#codegen
   %SEBFLUX using the new top node temperature, compute a new surface flux
   %
   % 
   % ea  - atmospheric vapor pressure from relative humidity data.
   % S   - stability function.
   % es  - surface saturation water vapor pressure.
   % Qe  - turbulent latent heat flux.
   % Qh  - turbulent sensible heat flux.
   % Qle - surface emitted longwave heat flux.
   % Qc  - conductive heat flux into the surface.
   % Qm  - the energy flux available for melting or freezing.
   % bal - an energy balance check.
   Tsfc     =  MELTTEMP(xTsfc,Tf);
   S        =  STABLEFN(Tair,Tsfc,wspd,scoef);
   es       =  VAPPRESS(Tsfc,Tf,liqflag);
   Qe       =  LATENT(De,S,ea,es,roL,epsilon,Pa);
   Qh       =  SENSIBLE(De,S,Tair,Tsfc,cv_air);
   Qle      =  LONGOUT(Tsfc,emiss,SB);
   Qc       =  CONDUCT(k_eff,T,dz,Tsfc);
   [Qm,Qf]  =  MFENERGY(albedo,Qsi,Qli,Qle,Qh,Qe,Qc,xTsfc,Tf,     ...
               Tair,wspd,De,ea,roL,Pa,cv_air,emiss,SB,k_eff,      ...
               T,dz,epsilon,scoef,chi);
   balance  =  ENBAL(albedo,Qsi,Qli,Qle,Qh,Qe,Qc,Qm,emiss,chi);
end
