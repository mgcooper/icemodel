function [Fsfc, Fdot] = SFCFLUX(Ta, Qsi, Qli, albedo, wspd, ppt, tppt, Pa, ...
      De, Ts, Qc, ea, cv_air, cv_liq, emiss, SB, roL, Tf, scoef, chi, liqflag)
   %SFCFLUX Evaluate the explicit bulk-Richardson surface residual.
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
   %
   %#codegen

   % Total surface flux
   Fsfc = surface_flux_residual(Ts, Ta, Qsi, Qli, albedo, wspd, ppt, tppt, ...
      Pa, De, Qc, ea, cv_air, cv_liq, emiss, SB, roL, scoef, chi, liqflag);
   
   % Numerical derivative
   Fdot = (surface_flux_residual(Ts + 1e-10, Ta, Qsi, Qli, albedo, wspd, ...
      ppt, tppt, Pa, De, Qc, ea, cv_air, cv_liq, emiss, SB, roL, scoef, ...
      chi, liqflag) - Fsfc) / 1e-10;

   % % for testing
   % Qr = EEE - emiss*SB.*Ts.^4;
   % Qh = AAA.*Sfnc(Ta,Ts,wspd,scoef).*(Ta-Ts)
   % Qe = FFF.*CCC.*Sfnc(Ta,Ts,wspd,scoef).*(ea-Vfnc(Ts,Tf,liqflag))
   %
   % % total heat = net radiation + sensible + latent
   % Q = Qr + Qh + Qe
end

function Fsfc = surface_flux_residual(Ts, Ta, Qsi, Qli, albedo, wspd, ppt, ...
      tppt, Pa, De, Qc, ea, cv_air, cv_liq, emiss, SB, roL, scoef, chi, ...
      liqflag)
   %SURFACE_FLUX_RESIDUAL Evaluate the explicit SEB residual used by SFCFLUX.
   %
   % This helper is file-local rather than nested because the residual has a
   % stable standalone contract and includes every flux term used by the
   % explicit bulk-Richardson surface solve, including precipitation advection.

   [Qe, Qh] = icemodel.surface.turbulent_heat_flux_bulk_richardson( ...
      Ta, Ts, wspd, Pa, De, ea, cv_air, roL, scoef, liqflag, NaN);
   Qa = QADVECT(ppt, tppt, cv_liq);

   Fsfc = chi * Qsi * (1.0 - albedo) + emiss * (Qli - SB * Ts ^ 4) ...
      + Qc + Qa + Qh + Qe;
end
