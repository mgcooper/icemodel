function [Fc, Fp] = surface_flux_linearization(T_sfc, tair, Qsi, Qli, ...
      albedo, wspd, ppt, tppt, ea_atm, H_h, H_e, br_coefs, liqflag, chi)
   %SURFACE_FLUX_LINEARIZATION Linearize the surface energy balance equation.
   %
   % The linearization is of the form: F = Fc + Fp * T
   %
   % To construct the linearization, the SEB is separated into terms that are
   % independent of T (call them F0) and terms that depend on T. The latter
   % are linearized using truncated Taylor expansions. Terms which are
   % independent of the current T (i.e., those that depend on past values of
   % T or the derivative dF/dT) are grouped with the F0 terms to form Fc:
   %
   %  SEB = F0 + F(T)
   %      = F0 + F' + (dF/dT)' * (T - T')
   %      = F0 + F' + (dF/dT)' * T - (dF/dT)' * T'
   %      = Fc + Fp * T
   %
   % where
   %
   %  Fc = F0 + F' - (dF/dT)' * T'
   %  Fp = (dF/dT)'
   %
   % Note that the outgoing longwave and the saturation vapor pressure are
   % linearized around T_old but the stability function is not linearized,
   % thus T_old should be used to compute the stability function.
   %
   % All terms passed here are 'old', meaning the linearizations are computed
   % once at the start of the timestep, and on iterations updated as
   % F = Fc + Fp * T_new
   %
   % Fc and Fp represent only the non-conductive surface flux linearization.
   % Conduction enters the Robin boundary-condition system through the top-node
   % finite-difference equation in `icemodel.column.assemble_enthalpy_system`
   % via the conductive conductance a1 = k_eff(1)/(dz(1)/2). In contrast, the
   % Dirichlet solve (solve_surface_temperature) includes dQc/dT_sfc directly in
   % the Newton-Raphson Jacobian.
   %
   % H_h  — sensible heat transport prefactor [W m-2 K-1] = cv_atm * De
   % H_e  — latent heat transport prefactor [W m-2 Pa-1] = hv_atm * De_e
   %
   % See also: icemodel.surface.numerical_surface_flux
   %
   %#codegen

   % Pass cv_liq to advective_heat_flux. TODO: For snowfall, pass snow density.
   persistent cv_liq
   if isempty(cv_liq)
      cv_liq = icemodel.physicalConstant('cv_liq');
   end

   % Surface saturation vapor pressure and derivative.
   [es_sfc, des_sfc_dT] = icemodel.vapor.saturation_vapor_pressure( ...
      T_sfc, liqflag);

   % Bulk richardson stability function.
   stability = ...
      icemodel.surface.turbulence.bulk_richardson.stability_factor( ...
      T_sfc, tair, wspd, br_coefs);

   % Turbulent heat fluxes via leaf functions (no diagnostics needed).
   Qe = ...
      icemodel.surface.turbulence.bulk_richardson.latent_heat_flux( ...
      es_sfc, ea_atm, H_e, stability);
   Qh = ...
      icemodel.surface.turbulence.bulk_richardson.sensible_heat_flux( ...
      T_sfc, tair, H_h, stability);

   % Canonical non-conductive surface terms at the current linearization
   % state. The Robin solve handles conduction in the interior column system.
   Qa = icemodel.surface.advective_heat_flux(ppt, tppt, cv_liq);
   Qsn = icemodel.surface.net_shortwave_radiation(Qsi, albedo, chi);
   [Qln, dQln_dT] = icemodel.surface.net_longwave_radiation(T_sfc, Qli);
   Q_sfc = icemodel.surface.evaluate_surface_energy_balance( ...
      Qsn, Qln, Qh, Qe, 0.0, Qa, 0.0);

   %%% Linearizations

   % Net longwave radiation.
   Fp_Qln = dQln_dT;

   % Sensible heat flux: H_h already absorbs cv_atm = ro_atm * cp_air
   Fp_Qh = -H_h * stability;

   % Latent heat flux (linearization of es_sfc around T_sfc):
   % H_e already absorbs hv_atm * De * epsilon / psfc
   Fp_Qe = -H_e * stability * des_sfc_dT;

   % Total linearization.
   Fp = Fp_Qln + Fp_Qh + Fp_Qe;
   Fc = Q_sfc - Fp * T_sfc;
end
