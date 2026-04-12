function [Q_sfc, dQ_sfc_dTs] = evaluate_surface_flux(tair, Qsi, Qli, ...
      albedo, wspd, ppt, tppt, psfc, De, T_sfc, ea_atm, roL, br_coefs, ...
      chi, liqflag, k_eff, T_ice, dz, ro_sfc, snow_depth, opts)
   %EVALUATE_SURFACE_FLUX Evaluate the bulk-Richardson SEB residual and derivative.
   %
   %  [Q_sfc, dQ_sfc_dTs] = ...
   %     icemodel.surface.turbulence.bulk_richardson.evaluate_surface_flux(...)
   %
   % This helper is a numerical-reference wrapper around the canonical
   % surface-energy-balance residual. It evaluates the bulk-Richardson
   % surface residual at `T_sfc` and estimates the derivative with a
   % complex-step perturbation:
   %
   %   dQ_sfc_dTs ≈ imag(Q_sfc(T_sfc + 1i*h)) / h,  h = 1e-10
   %
   % The helper exists primarily for derivative-validation tests; the
   % production Dirichlet solve still uses the analytical Jacobian in
   % `icemodel.surface.solve_surface_temperature`.
   % The passed opts struct is copied locally with
   % `turbulent_flux_scheme = 'bulk_richardson'` so this wrapper always
   % evaluates the bulk-Richardson residual regardless of the caller's
   % active scheme setting.
   %
   % See also: icemodel.surface.surface_energy_balance_residual,
   %           icemodel.surface.solve_surface_temperature,
   %           icemodel.numerics.complex_step_derivative
   %
   %#codegen

   opts_eval = opts;
   opts_eval.turbulent_flux_scheme = 'bulk_richardson';

   residual_fn = @(Ts_eval) icemodel.surface.surface_energy_balance_residual( ...
      Ts_eval, tair, Qsi, Qli, albedo, wspd, ppt, tppt, psfc, De, ea_atm, ...
      br_coefs, roL, liqflag, chi, T_ice, k_eff, dz, ro_sfc, snow_depth, opts_eval);

   Q_sfc = residual_fn(T_sfc);
   dQ_sfc_dTs = icemodel.numerics.complex_step_derivative(residual_fn, T_sfc);
end
