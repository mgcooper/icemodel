function [Fsfc, Fdot] = surface_fluxes(tair, Qsi, Qli, albedo, wspd, ppt, ...
      tppt, psfc, De, T_sfc, Qc, ea_atm, roL, br_coefs, chi, liqflag)
   %SURFACE_FLUXES Evaluate the explicit bulk-Richardson surface residual.
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
   Fsfc = surface_flux_residual(T_sfc, tair, Qsi, Qli, albedo, wspd, ppt, tppt, ...
      psfc, De, Qc, ea_atm, roL, br_coefs, chi, liqflag);

   % Numerical derivative
   Fdot = (surface_flux_residual(T_sfc + 1e-10, tair, Qsi, Qli, albedo, wspd, ...
      ppt, tppt, psfc, De, Qc, ea_atm, roL, br_coefs, ...
      chi, liqflag) - Fsfc) / 1e-10;

   % % for testing
   % Qr = EEE - emiss*SB.*Ts.^4;
   % Qh = AAA.*Sfnc(Ta,Ts,wspd,scoef).*(Ta-Ts)
   % Qe = FFF.*CCC.*Sfnc(Ta,Ts,wspd,scoef).*(ea-Vfnc(Ts,Tf,liqflag))
   %
   % % total heat = net radiation + sensible + latent
   % Q = Qr + Qh + Qe
end

function Fsfc = surface_flux_residual(T_sfc, tair, Qsi, Qli, albedo, wspd, ppt, ...
      tppt, psfc, De, Qc, ea_atm, roL, br_coefs, chi, liqflag)
   %SURFACE_FLUX_RESIDUAL Evaluate the explicit SEB residual used here.
   %
   % This helper is file-local rather than nested because the residual has a
   % stable standalone contract and includes every flux term used by the
   % explicit bulk-Richardson surface solve, including precipitation advection.

   persistent cv_liq emiss SB
   if isempty(cv_liq)
      [cv_liq, SB] = icemodel.physicalConstant('cv_liq', 'SB');
      emiss = icemodel.parameterLookup('emiss');
   end

   [Qe, Qh] = ...
      icemodel.surface.turbulence.bulk_richardson.turbulent_heat_flux( ...
      T_sfc, tair, wspd, psfc, ea_atm, De, br_coefs, roL, liqflag, NaN);
   Qa = QADVECT(ppt, tppt, cv_liq);

   Fsfc = chi * Qsi * (1.0 - albedo) + emiss * (Qli - SB * T_sfc ^ 4) ...
      + Qc + Qa + Qh + Qe;
end
