function [Qe, Qh, Qc, Qm, Qf, balance] = SEBFLUX(T, xTs, Ta, Qsi, Qli, ...
      albedo, wspd, ppt, tppt, Pa, De, ea, Tf, k_eff, dz, cv_air, cv_liq, ...
      roL, emiss, SB, chi, epsilon, scoef, liqflag)
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
   %
   %#codegen

   Ts       =  MELTTEMP(xTs,Tf);
   S        =  STABLEFN(Ta,Ts,wspd,scoef);
   es       =  VAPPRESS(Ts,Tf,liqflag);
   Qe       =  LATENT(De,S,ea,es,roL,epsilon,Pa);
   Qh       =  SENSIBLE(De,S,Ta,Ts,cv_air);
   Qle      =  LONGOUT(Ts,emiss,SB);
   Qc       =  CONDUCT(k_eff,T,dz,Ts);
   Qa       =  QADVECT(ppt,tppt,cv_liq);
   [Qm,Qf]  =  MFENERGY(albedo,Qsi,Qli,Qle,Qh,Qe,Qc,Qa,xTs,Tf, ... % Note: xTs
      Ta,wspd,ppt,tppt,De,ea,roL,Pa,cv_liq,cv_air,emiss,SB,k_eff, ...
      T,dz,epsilon,scoef,chi);
   balance  =  ENBAL(albedo,emiss,chi,Qsi,Qli,Qle,Qh,Qe,Qc,Qa,Qm);
end
