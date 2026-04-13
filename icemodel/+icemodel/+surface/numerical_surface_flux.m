function [Q_sfc, dQ_sfc_dTs] = numerical_surface_flux(T_sfc, tair, Qsi, ...
      Qli, albedo, wspd, ppt, tppt, psfc, De, ea_atm, br_coefs, ro_air_Lv, ...
      liqflag, chi, T_ice, k_eff, dz, ro_sfc, snow_depth, opts)
   %NUMERICAL_SURFACE_FLUX Evaluate the SEB residual and derivative numerically.
   %
   %  [Q_sfc, dQ_sfc_dTs] = icemodel.surface.numerical_surface_flux(...)
   %
   % This helper is the canonical numerical-reference wrapper around the
   % surface-energy-balance residual. It evaluates the current SEB residual
   % at `T_sfc` using the caller's active turbulent-flux scheme in `opts`,
   % then estimates the derivative with a complex-step perturbation:
   %
   %   dQ_sfc_dTs ≈ imag(Q_sfc(T_sfc + 1i*h)) / h,  h = 1e-10
   %
   % The helper exists primarily for derivative-validation tests; the
   % production Dirichlet solve still uses the analytical Jacobian in
   % `icemodel.surface.solve_surface_temperature`, and the production Robin
   % linearization may still use a scheme-specific analytical helper.
   %
   % See also: icemodel.surface.surface_energy_balance_residual,
   %           icemodel.surface.surface_flux_linearization,
   %           icemodel.surface.solve_surface_temperature,
   %           icemodel.numerics.complex_step_derivative
   %
   %#codegen

   residual_fn = @(Ts_eval) ...
      icemodel.surface.surface_energy_balance_residual(Ts_eval, tair, Qsi, ...
      Qli, albedo, wspd, ppt, tppt, psfc, De, ea_atm, br_coefs, ro_air_Lv, ...
      liqflag, chi, T_ice, k_eff, dz, ro_sfc, snow_depth, opts);

   Q_sfc = residual_fn(T_sfc);
   dQ_sfc_dTs = icemodel.numerics.complex_step_derivative(residual_fn, T_sfc);
end
