function residual = surface_energy_balance_residual(T_sfc, tair, Qsi, Qli, ...
      albedo, wspd, ppt, tppt, psfc, De, ea_atm, br_coefs, ro_air_Lv, liqflag, ...
      chi, T_ice, k_eff, dz, ro_sfc, snow_depth, opts)
   %SURFACE_ENERGY_BALANCE_RESIDUAL Return the SEB residual at T_sfc.
   %
   %  residual = icemodel.surface.surface_energy_balance_residual(...)
   %
   % This is the canonical nonlinear surface residual used by the touched
   % SEB solvers. It wraps the surface term diagnostics and then evaluates
   % the residual with no melt term.
   %
   %#codegen

   [Qe, Qh, Qc, Qa, Qsn, Qln] = icemodel.surface.surface_energy_balance_terms( ...
      T_sfc, tair, Qsi, Qli, albedo, wspd, ppt, tppt, psfc, De, ea_atm, ...
      br_coefs, ro_air_Lv, liqflag, chi, T_ice, k_eff, dz, ro_sfc, snow_depth, ...
      opts);

   residual = icemodel.surface.evaluate_surface_energy_balance( ...
      Qsn, Qln, Qh, Qe, Qc, Qa, 0.0);
end
