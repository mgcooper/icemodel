%------------------------------------------------------------------------------
%   COMPUTE THE SURFACE ENERGY-BALANCE AND SOLVE FOR Tsfc
%------------------------------------------------------------------------------
function [  Qm,                                                         ...
            Qf,                                                         ...
            Qh,                                                         ...
            Qe,                                                         ...
            Qc,                                                         ...
            Qle,                                                        ...
            balance,                                                    ...
            Tsfc  ]  =  ENBALANCE(Tair,wspd,rh,Qsi,Qli,albedo,k_eff,Pa, ...
                        T,Tf,dz,chi,xTsfc,cv_air,emiss,SB,roL,De,scoef, ...
                        epsilon,fopts,liqflag,opts)
%------------------------------------------------------------------------------
   
% Atmospheric vapor pressure from relative humidity data.
   ea       =  VAPPRESS(rh,Tair,liqflag);
   
% incoming longwave if not provided
%  Qli      =  LONGIN(ea,Tair,stefBoltz);
   
% Compute the average station pressure.
%  Pa       =  PRESSURE(topo);
   
% Compute the turbulent exchange coefficients.
%  De       =  EXCOEFS(wspd,wcoef);
   
% Compute the flux contribution due to conduction.
   Qc       =  CONDUCT(k_eff,T,dz,xTsfc);
   
% Solve the energy balance for the surface temperature.
   [Tsfc,~] =  SFCTEMP(Tair,Qsi,Qli,ea,albedo,De,Pa,wspd,cv_air,      ...
               emiss,SB,Tf,Qc,xTsfc,chi,roL,scoef,fopts,liqflag);
   
% Make the Tsfc_0 <= 0 C for surface flux calculations.
%   Let Tsfc remain > Tf for the upper boundary condition on ICE_ENERGY
%    Tsfc0    =  MELTTEMP(Tsfc,Tf);
   
% Compute the stability function.
   S        =  STABLEFN(Tair,MELTTEMP(Tsfc,Tf),wspd,scoef);
   
% Compute the water vapor pressure at the surface.
   es0      =  VAPOR(MELTTEMP(Tsfc,Tf),Tf,liqflag);
   
% Compute the latent heat flux.
   Qe       =  LATENT(De,S,ea,es0,roL,epsilon,Pa);

% Compute the sensible heat flux.
   Qh       =  SENSIBLE(De,S,Tair,MELTTEMP(Tsfc,Tf),cv_air);

% Compute the longwave flux emitted by the surface.
   Qle      =  LONGOUT(MELTTEMP(Tsfc,Tf),emiss,SB);
   
% Compute the energy flux available for melting or freezing.
   Qm       =  0.0; 
   Qf       =  0.0;

% if melting, compute melt energy   
   if Tsfc>=Tf 
      Qm    =  chi*Qsi*(1.0-albedo) + emiss*Qli + Qle + Qh + Qe + Qc;
   else
% else compute energy needed to reach melt temp   
      Qf    =  -(chi*(1.0-albedo)*Qsi+emiss*Qli+LONGOUT(Tf,emiss,SB)+...
               LATENT(De,STABLEFN(Tair,Tf,wspd,scoef),ea,...
               VAPOR(Tf,Tf,true),roL,epsilon,Pa)+...
               SENSIBLE(De,STABLEFN(Tair,Tf,wspd,scoef),Tair,Tf,cv_air)+...
               CONDUCT(k_eff,T,dz,Tf));
   end

% Perform an energy balance check.
   balance  =  chi*Qsi*(1.0-albedo) + emiss*Qli + Qle + Qh + Qe + Qc - Qm;
   
% For a 'skin' surface energy balance model, reset Tsfc
   if strcmp('skinmodel', opts.simmodel)
      Tsfc  =  Tsfc0;
   end
