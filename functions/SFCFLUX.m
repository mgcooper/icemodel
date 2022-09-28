%--------------------------------------------------------------------------
%   SOLVE THE ENERGY BALANCE FOR SURFACE TEMPERATURE
%--------------------------------------------------------------------------

function [Fsfc,Fdot] =  SFCFLUX(Tair,Qli,Qc,ea,De,Pa,wspd,cv_air,emiss, ...
                        SB,roL,Tf,Tsfc,scoef,liqflag,chi,Qsi,albedo)
%--------------------------------------------------------------------------

% The surface flux includes Qc so it can be used in general, but for a
% neumann upper bc at the top, Qc=0 when this is
% called from within the solver. Tsfc is then updated using the top node
% temperature and the conductive flux from the top node into the surface.
% In contrast, for a dirichlet upper bc, this function could be called with
% a known Qc (using the previous timestep for instance) to solve for Tsfc,
% which would then be fed into the solver. 

% regarding the source term, if chi=0, then no solar heat is included in
% FSFC except the portion that is already included in Sc calculated in
% SOLAR_HEAT. if chi~=0, then the portion allocated to the 'skin' is
% included. This portion represents the longwave energy that does not
% penetrate the surface more than a mm at most. 

% gather terms in the SEB equation. Note that Qc = 0 for Nuemann bc.
   AAA   =  cv_air * De;                           % [W m-2 K-1]
   CCC   =  0.622 / Pa;                            % [Pa-1] = [m3 J-1]
   EEE   =  chi*Qsi*(1.0-albedo) + emiss*Qli + Qc;     % [W m-2]
   FFF   =  roL * De;                             % [W m-2]
   
% % Compute the constants used in the stability coefficient computations
%    B1    =  scoef(2)/(Tair*wspd^2);                % [K-1]
%    B2    =  scoef(3)/(sqrt(Tair)*wspd);            % [K-1]

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% fzero
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   Sfnc  =  @STABLEFN;
   Vfnc  =  @VAPOR;
   fSEB  =  @(Tsfc) EEE - emiss*SB.*Tsfc.^4 +                            ...
               AAA.*Sfnc(Tair,Tsfc,wspd,scoef).*(Tair-Tsfc)             ...
             + FFF.*CCC.*Sfnc(Tair,Tsfc,wspd,scoef)                     ...
            .* (ea-Vfnc(Tsfc,Tf,liqflag));
             %  + cp_liq*ppt*Tppt;                % ppt in kg/m2/s
   Fsfc    =   fSEB(Tsfc);
   Fdot    =   (fSEB(Tsfc+1e-10)-Fsfc)/1e-10;
    
% % for testing    
% Qr =   EEE - emiss*SB.*Tsfc.^4;
% Qh =   AAA.*Sfnc(Tair,Tsfc,wspd,scoef).*(Tair-Tsfc)
% Qe =   FFF.*CCC.*Sfnc(Tair,Tsfc,wspd,scoef).*(ea-Vfnc(Tsfc,Tf,liqflag))
%     
% % total heat = net radiation + sensible + latent
% Q = Qr + Qh + Qe
    
    
    