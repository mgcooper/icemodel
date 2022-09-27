function [  Qe,                                                         ...
            Qh,                                                         ...
            Qc,                                                         ...
            Qm,                                                         ...
            Qf,                                                         ...
            ea,                                                         ...
            balance] = SEBFLUX(T,Tair,Tsfc,Qsi,Qli,albedo,wspd,rh,Pa,   ...
                       Tf,k_eff,cv_air,roL,emiss,SB,epsilon,scoef,De,   ...
                       dz,liqflag,chi)

% using the new top node temperature, compute a new surface flux
                        
% Atmospheric vapor pressure from relative humidity data.
% Compute the stability function.   
% Compute the water vapor pressure at the surface.
% Compute the latent heat flux.
% Compute the sensible heat flux.
% Compute the longwave flux emitted by the surface.
% Compute the flux contribution due to conduction.
% Compute the energy flux available for melting or freezing.
% Perform an energy balance check.
   ea       =  VAPPRESS(rh,Tair,liqflag);
   S        =  STABLEFN(Tair,Tsfc,wspd,scoef);
   es0      =  VAPOR(Tsfc,Tf,liqflag);
   Qe       =  LATENT(De,S,ea,es0,roL,epsilon,Pa);
   Qh       =  SENSIBLE(De,S,Tair,Tsfc,cv_air);
   Qle      =  LONGOUT(Tsfc,emiss,SB);
   Qc       =  CONDUCT(k_eff,T,dz,Tsfc);
   [Qm,Qf]  =  MFENERGY(albedo,Qsi,Qli,Qle,Qh,Qe,Qc,Tsfc,Tf,      ...
               Tair,wspd,De,ea,roL,Pa,cv_air,emiss,SB,k_eff,      ...
               T,dz,epsilon,scoef,chi);
   balance  =  ENBAL(albedo,Qsi,Qli,Qle,Qh,Qe,Qc,Qm,emiss,chi);
   
   
   