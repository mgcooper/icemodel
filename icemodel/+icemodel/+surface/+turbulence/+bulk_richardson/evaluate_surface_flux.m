function [Q_sfc, dQ_sfc_dTs] = evaluate_surface_flux(tair, Qsi, Qli, ...
      albedo, wspd, ppt, tppt, psfc, De, T_sfc, Qc, ea_atm, roL, ...
      br_coefs, chi, liqflag)
   %EVALUATE_SURFACE_FLUX Evaluate the bulk-Richardson SEB residual and derivative.
   %
   %  [Q_sfc, dQ_sfc_dTs] = ...
   %     icemodel.surface.turbulence.bulk_richardson.evaluate_surface_flux(...)
   %
   % Evaluates the explicit surface energy balance at T_sfc:
   %
   %   Q_sfc = chi*Qsi*(1-albedo) + emiss*(Qli - SB*T_sfc^4)
   %           + Qh(T_sfc) + Qe(T_sfc) + Qc + Qa
   %
   % and returns its numerical derivative using a Newton difference quotient:
   %
   %   dQ_sfc_dTs ≈ (Q_sfc(T_sfc + h) - Q_sfc(T_sfc)) / h,  h = 1e-10
   %
   % --- Relation to solve_surface_temperature ---
   %
   % icemodel.surface.solve_surface_temperature minimizes Q_sfc via
   % Newton-Raphson using an analytical Jacobian. evaluate_surface_flux
   % provides an independent numerical evaluation of the same residual and
   % its derivative, which serves as:
   %   (1) A validation tool: verify that solve_surface_temperature converged
   %       to a state where Q_sfc is small (see test_sfctemp_finds_small_residual).
   %   (2) A derivative check: verify the analytical Jacobian in
   %       solve_surface_temperature agrees with the numerical one
   %       (see test_sfcflux_derivative_matches_finite_difference).
   %
   % The analytical Jacobian is kept in solve_surface_temperature rather than
   % delegating to this function's dQ_sfc_dTs for two reasons:
   %   (a) Efficiency: the analytical form needs one evaluation per Newton step,
   %       while the numerical form needs two (F(T) and F(T+h)).
   %   (b) Exactness: the analytical form includes the conduction Jacobian term
   %       -a1 = -k_eff(1)/(dz(1)/2) when Qc couples back to T_sfc; the
   %       numerical form here treats Qc as fixed and therefore misses this
   %       term when Qc is column-state-dependent.
   %
   % --- Relation to surface_flux_linearization ---
   %
   % bulk_richardson/surface_flux_linearization returns [Fc, Fp] where
   % Q_sfc ≈ Fc + Fp*T_sfc, using an analytical Taylor expansion. It is used
   % for the Robin boundary condition coupling and intentionally excludes Qc
   % (which is handled by the interior subsurface solve).
   %
   % evaluate_surface_flux includes Qc and returns [Q_sfc, dQ_sfc_dTs] —
   % the residual and its raw numerical derivative rather than the linearized
   % form. The two functions serve distinct solver paths.
   %
   % See also: icemodel.surface.solve_surface_temperature,
   %           icemodel.surface.turbulence.bulk_richardson.surface_flux_linearization,
   %           icemodel.surface.surface_energy_balance_residual
   %
   %#codegen

   persistent h
   if isempty(h)
      h = 1e-10;
   end

   Q_sfc = surface_flux(T_sfc, tair, Qsi, Qli, albedo, wspd, ppt, tppt, ...
      psfc, De, Qc, ea_atm, roL, br_coefs, chi, liqflag);

   dQ_sfc_dTs = (surface_flux(T_sfc + h, tair, Qsi, Qli, albedo, wspd, ppt, ...
      tppt, psfc, De, Qc, ea_atm, roL, br_coefs, chi, liqflag) - Q_sfc) / h;
end

function Q_sfc = surface_flux(T_sfc, tair, Qsi, Qli, albedo, wspd, ppt, ...
      tppt, psfc, De, Qc, ea_atm, roL, br_coefs, chi, liqflag)
   %SURFACE_FLUX Evaluate the non-conductive + conductive SEB residual at T_sfc.
   %
   % This is the explicit surface flux evaluated at a fixed Qc. Used by the
   % public evaluate_surface_flux at two temperatures to form the difference
   % quotient.

   persistent cv_liq emiss SB
   if isempty(cv_liq)
      [cv_liq, SB] = icemodel.physicalConstant('cv_liq', 'SB');
      emiss = icemodel.parameterLookup('emiss');
   end

   [Qe, Qh] = ...
      icemodel.surface.turbulence.bulk_richardson.turbulent_heat_flux( ...
      T_sfc, tair, wspd, psfc, ea_atm, De, br_coefs, roL, liqflag, NaN);
   Qa = icemodel.surface.advective_heat_flux(ppt, tppt, cv_liq);

   Q_sfc = chi * Qsi * (1.0 - albedo) + emiss * (Qli - SB * T_sfc ^ 4) ...
      + Qc + Qa + Qh + Qe;
end
