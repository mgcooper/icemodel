function [T_sfc, ok] = solve_surface_temperature(T_sfc, tair, Qsi, Qli, ...
      albedo, wspd, ppt, tppt, ea_atm, H_h, H_e, br_coefs, liqflag, ...
      chi, T_ice, k_eff, dz)
   %SOLVE_SURFACE_TEMPERATURE Solve the explicit bulk-Richardson SEB for Ts.
   %
   %  [T_sfc, ok] = icemodel.surface.solve_surface_temperature(T_sfc, tair, ...)
   %
   % Solves the nonlinear surface energy balance for T_sfc using a
   % Newton-Raphson iteration with an analytical Jacobian. The residual
   % being minimized is:
   %
   %   F(T_sfc) = chi*Qsi*(1-albedo) + emiss*(Qli - SB*T_sfc^4)
   %              + Qh(T_sfc) + Qe(T_sfc) + Qc(T_sfc) + Qa
   %
   % where Qc(T_sfc) = k_eff(1) * (T_ice(1) - T_sfc) / (dz(1)/2) couples
   % the surface temperature to the top ice layer, and Qh, Qe are evaluated
   % using the bulk-Richardson turbulence closure. The analytical Jacobian is:
   %
   %   dF/dT_sfc = -4 * emiss * SB * T_sfc^3
   %               + d(stability)/dT * (Qh/stability + Qe/stability)
   %               + stability * (dQh/dT|stability + dQe/dT|stability)
   %               + dQc/dT_sfc
   %
   % where dQc/dT_sfc = -k_eff(1) / (dz(1) / 2) is constant.
   %
   % H_h and H_e are the precomputed sensible and latent heat transport
   % prefactors from initialize_surface_state / updatesubstep.
   %
   % This is the Dirichlet surface solve: Qc and its derivative enter the
   % Newton-Raphson residual and Jacobian directly. In the Robin path,
   % conduction instead enters through the top-node finite-difference
   % equation in `icemodel.column.assemble_enthalpy_system` rather than as
   % an explicit Jacobian term here.
   %
   % T_sfc on input is the outer coupling iterate from
   % solve_surface_energy_balance, used as the initial Newton guess. On
   % convergence T_sfc is the solution.
   %
   % See also: icemodel.surface.solve_surface_energy_balance,
   %           icemodel.surface.conductive_heat_flux,
   %           icemodel.surface.surface_energy_balance_residual
   %
   %#codegen

   % Solver tolerances.
   persistent tol maxiter
   if isempty(tol)
      tol = 1e-3;
      maxiter = 100;
   end

   % Pass cv_liq to advective_heat_flux. TODO: For snowfall, pass snow density.
   persistent cv_liq
   if isempty(cv_liq)
      cv_liq = icemodel.physicalConstant('cv_liq');
   end

   %%% Gather T_sfc-independent terms in the SEB equation.
   Qsn = icemodel.surface.net_shortwave_radiation(Qsi, albedo, chi);
   Qa = icemodel.surface.advective_heat_flux(ppt, tppt, cv_liq);

   % Derivative of conductive heat flux is T_sfc-independent.
   [~, dQc_dT_sfc] = icemodel.surface.conductive_heat_flux( ...
      k_eff, T_ice, dz, T_sfc);

   ok = false;
   old = T_sfc;
   for iter = 1:maxiter

      %%% Compute heat fluxes and Jacobian terms at the current T_sfc iterate.

      % Surface saturation vapor pressure and its temperature derivative.
      [es_sfc, des_sfc_dT] = icemodel.vapor.saturation_vapor_pressure( ...
         old, liqflag);

      % Bulk-Richardson stability factor and its temperature derivative.
      [stability, dstability] = ...
         icemodel.surface.turbulence.bulk_richardson.stability_factor( ...
         old, tair, wspd, br_coefs);

      % Latent heat flux and its temperature derivative.
      [Qe, dQe_dT_sfc] = ...
         icemodel.surface.turbulence.bulk_richardson.latent_heat_flux( ...
         es_sfc, ea_atm, H_e, stability, des_sfc_dT, dstability);

      % Sensible heat flux and its temperature derivative.
      [Qh, dQh_dT_sfc] = ...
         icemodel.surface.turbulence.bulk_richardson.sensible_heat_flux( ...
         old, tair, H_h, stability, dstability);

      % Net longwave radiation and its temperature derivative.
      [Qln, dQln_dTsfc] = icemodel.surface.net_longwave_radiation(old, Qli);

      % Conductive heat flux.
      Qc = icemodel.surface.conductive_heat_flux(k_eff, T_ice, dz, old);

      % Newton-Raphson residual and analytical Jacobian.
      f = icemodel.surface.evaluate_surface_energy_balance( ...
         Qsn, Qln, Qh, Qe, Qc, Qa, 0.0);

      dfdT = dQln_dTsfc + dQh_dT_sfc + dQe_dT_sfc + dQc_dT_sfc;

      % Updated T_sfc iterate.
      T_sfc = old - f / dfdT;

      if abs(T_sfc - old) < tol
         ok = true;
         return
      elseif ~isfinite(T_sfc) || iter == maxiter
         ok = false;
         T_sfc = tair;
         return
      end
      old = T_sfc;
   end
end
