function [Fsfc, Fdot] = SFCFLUX(Ta, Qsi, Qli, albedo, wspd, Pa, De, ...
      Ts, Qc, ea, cv_air, emiss, SB, roL, Tf, scoef, chi, liqflag)
   %SFCFLUX Solve the energy balance for surface temperature
   %
   % The surface flux includes Qc for general use, but note for a
   % coupled surface/column with neumann bc at the top, Qc = 0 when this is
   % called from within the solver. Ts is then updated using the top node
   % temperature and the conductive flux from the top node into the surface.
   % In contrast, for a dirichlet upper bc, this function could be called with
   % a known Qc (using the previous timestep for instance) to solve for Ts,
   % which would then be fed into the solver.
   %
   % Regarding the source term, if chi = 0, then no solar heat is included in
   % FSFC except the portion that is already included in Sc calculated in
   % SOLAR_HEAT. if chi ~= 0, then the portion allocated to the 'skin' is
   % included. This portion represents the longwave energy that does not
   % penetrate the surface more than a mm at most.

   % gather terms in the SEB equation. Note that Qc = 0 for Nuemann bc.
   AAA = cv_air * De;                              % [W m-2 K-1]
   CCC = 0.622 / Pa;                               % [Pa-1] = [m3 J-1]
   EEE = chi*Qsi*(1.0-albedo) + emiss*Qli + Qc;    % [W m-2]
   FFF = roL * De;                                 % [W m-2]

   % % Compute the constants used in the stability coefficient computations
   % B1 = scoef(2)/(Ta*wspd^2);                  % [K-1]
   % B2 = scoef(3)/(sqrt(Ta)*wspd);              % [K-1]

   % fzero
   Sfnc = @STABLEFN;
   Vfnc = @VAPPRESS;
   fSEB = @(Ts) ...
      EEE - emiss*SB.*Ts.^4 + ...
      AAA.*Sfnc(Ta,Ts,wspd,scoef).*(Ta-Ts) ...
      + FFF.*CCC.*Sfnc(Ta,Ts,wspd,scoef) ...
      .* (ea-Vfnc(Ts,Tf,liqflag));
   % + cp_liq*ppt*Tppt; % ppt in kg/m2/s
   Fsfc = fSEB(Ts);
   Fdot = (fSEB(Ts+1e-10)-Fsfc)/1e-10;

   % % for testing
   % Qr = EEE - emiss*SB.*Ts.^4;
   % Qh = AAA.*Sfnc(Ta,Ts,wspd,scoef).*(Ta-Ts)
   % Qe = FFF.*CCC.*Sfnc(Ta,Ts,wspd,scoef).*(ea-Vfnc(Ts,Tf,liqflag))
   %
   % % total heat = net radiation + sensible + latent
   % Q = Qr + Qh + Qe
end
