%--------------------------------------------------------------------------
% linearize the surface energy balance equation in the form S = Sc + Sp*T
%--------------------------------------------------------------------------

function [Sc,Sp]  =  SFCFLIN2(Tair,Qli,ea,De,Pa,wspd,cv_air,emiss,SB,    ...
                     roL,Tf,Tsfc,scoef,liqflag,chi,Qsi,albedo)
%--------------------------------------------------------------------------

% this is a different linearization that linearizes the outgoing longwave
% and the saturation vapor pressure around T(1)_old but does not linearize
% the stability function therefore T(1)_old has to be used to compute the
% stability function

% NOTE: all terms here are 'old', meaning these are computed once at the
% start of the timestep, and on iterations updated as S = Sc + Sp*T_new

% the linearization is of the form: SEB=SEB'+dSEB/dTsfc(Tsfc-Tsfc')

% note on solar radiation source term. the penetrating solar radiation
% absorbed in each cv (dQp) that comes out of updateextcoefs is scaled to
% match the total subsurface absorbed solar radiation at each timestep in
% SOLAR_HEAT (dQp_new = -(1-chi)*Qsi/total_solar).*dQp_old). Although not
% intuitive, just trust that this provides the amount of solar radiation
% absorbed in each cv AFTER removing chi*Qsi*1-albedo from Qsi. TLDR: as
% long as the SEB term gets chi*(1.0-albedo) * Qsi, and the Qsip that comes
% out of SOLARPEN and goes into SOLAR_HEAT uses the same chi, then energy
% is not double counted. The 'chi' part goes into the skin, and the rest
% goes into each cv. So, in the SEB linearization, we don't need any
% special procedures to deal with the part absorbed in the top layer. The
% equation below puts the skin part into the SEB linearization, and the
% heat solver puts in the cv part, and they are combined in the equation
% for for the top layer/node. 

if liqflag == true
   % Over water.
   A  =  611.210;
   B  =  17.502;
   C  =  240.97;
else
   % Over ice.
   A  =  611.150;
   B  =  22.452;
   C  =  272.55;
end

% Compute the constants used in the stability coefficient computations
   B1    =  scoef(2)/(Tair*wspd^2);
   B2    =  scoef(3)/(sqrt(Tair)*wspd);
   
% This accounts for an increase in turbulent fluxes under unstable conditions.
   es0   =  A*exp((B*(Tsfc-Tf))/(C+(Tsfc-Tf)));
    
   if (Tsfc>Tair)         % Unstable case.
      B3    =   1 + B2*sqrt(Tsfc-Tair);
      S     =   1 + B1*(Tsfc-Tair)/B3;
   elseif (Tsfc<Tair)     % Stable case.
      S     =   1 /(1 + B1/2*(Tair-Tsfc))^2;
   else                    % Neutrally stable case.
      S     =   1.0;
   end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% linearizations
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Note that Qc = 0 for Nuemann bc.

% outgoing longwave:
   Sc_Qle   =  3*emiss*SB*Tsfc^4;
   Sp_Qle   =  -4/3*Sc_Qle/Tsfc;

% sensible heat flux:
%    Sc_Qh    =  cv_air*De*(Tsfc*dSdT*(Tsfc-Tair) + Tair*S);
%    Sp_Qh    =  cv_air*De*(dSdT*(Tair-Tsfc)-S);
   Sc_Qh    =  cv_air*De*S*Tair;
   Sp_Qh    =  -cv_air*De*S;
   
% latent heat flux:
%    Sc_Qe    =  roL*De*0.622/Pa*(Tsfc*(dSdT*(es0-ea)+S*es0*AS)+S*(ea-es0));
%    Sp_Qe    =  roL*De*0.622/Pa*(dSdT*(ea-es0)-S*es0*AS);
   Sc_Qe    =  roL*De*0.622/Pa*S*(ea-es0*(1-B*C*Tsfc/(C+Tsfc-Tf)^2));
   Sp_Qe    =  -roL*De*0.622/Pa*S*es0*B*C/(C+Tsfc-Tf)^2;
   
% combine net sw, incoming lw, conduction, and snow/rain heat flux:
   Sc       =  Sc_Qle + Sc_Qh + Sc_Qe + emiss*Qli + chi*Qsi*(1.0-albedo); 
   Sp       =  Sp_Qle + Sp_Qh + Sp_Qe;

% might need to add a check for switch from stable/unstable within a step
% check if Sp_Qh or Sp_Qe should be negative

    
% % for testing    
% Qr =   EEE - emiss*SB.*xTsfc.^4;
% Qh =   AAA.*Sfnc(Tair,xTsfc,wspd,grav,kapp,z_0,z_obs).*(Tair-xTsfc)
% Qe =   FFF.*CCC.*Sfnc(Tair,xTsfc,wspd,grav,kapp,z_0,z_obs).*(ea-Vfnc(xTsfc,Tf))
%     
% % total heat = net radiation + sensible + latent
% Q = Qr + Qh + Qe
    
% don't think this is complete might have been distracted before finishing
% SEB =  EEE+3*emiss*SB*xTsfc^4+Sc_Qh+Sc_Qe+(-4*emiss*SB*T_old(1)^3-Sp_Q
    
    