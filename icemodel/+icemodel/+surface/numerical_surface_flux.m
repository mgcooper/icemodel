function [Q_sfc, dQ_sfc_dTs] = numerical_surface_flux(T_sfc, tair, Qsi, ...
      Qli, albedo, wspd, ppt, tppt, psfc, ea_atm, ro_atm, cv_atm, nu_air, ...
      H_h, H_e, hv_atm, br_coefs, liqflag, chi, T_ice, k_eff, dz, ro_sfc, ...
      snow_depth, opts)
   %NUMERICAL_SURFACE_FLUX Evaluate the SEB residual and derivative numerically.
   %
   %  [Q_sfc, dQ_sfc_dTs] = icemodel.surface.numerical_surface_flux(...)
   %
   % Evaluates the current SEB residual at T_sfc using the turbulent-flux
   % scheme configured in opts, then estimates the derivative with a
   % complex-step perturbation:
   %
   %   dQ_sfc_dTs ≈ imag(Q_sfc(T_sfc + 1i*h)) / h,  h = 1e-10
   %
   % The function exists mainly for derivative-validation tests. The Dirichlet
   % SEB solve uses the analytical Jacobian in
   % icemodel.surface.solve_surface_temperature, and the Robin linearization may
   % use a scheme-specific analytical helper.
   %
   % See also: icemodel.surface.surface_energy_balance_residual,
   %           icemodel.surface.surface_flux_linearization,
   %           icemodel.surface.solve_surface_temperature,
   %           icemodel.numerics.complexstep_derivative
   %
   %#codegen

   residual_fn = @(Ts_eval) ...
      icemodel.surface.surface_energy_balance_residual(Ts_eval, tair, ...
      Qsi, Qli, albedo, wspd, ppt, tppt, psfc, ea_atm, ro_atm, ...
      cv_atm, nu_air, H_h, H_e, hv_atm, br_coefs, liqflag, chi, ...
      T_ice, k_eff, dz, ro_sfc, snow_depth, opts);

   Q_sfc = residual_fn(T_sfc);
   dQ_sfc_dTs = icemodel.numerics.complexstep_derivative( ...
      residual_fn, T_sfc);
end
