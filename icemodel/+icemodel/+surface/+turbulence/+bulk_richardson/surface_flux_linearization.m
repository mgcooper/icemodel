function [Fc, Fp] = surface_flux_linearization(T_sfc, tair, Qsi, Qli, albedo, ...
      wspd, ppt, tppt, psfc, De, ea_atm, br_coefs, ro_air_Lv, liqflag, chi)
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
   % See also: icemodel.surface.numerical_surface_flux
   %
   %#codegen

   persistent cv_air cv_liq emiss SB epsilon
   if isempty(cv_air)
      [cv_air, cv_liq, SB, epsilon] = icemodel.physicalConstant( ...
         'cv_air', 'cv_liq', 'SB', 'epsilon');
      emiss = icemodel.parameterLookup('emiss');
   end

   % Surface saturation vapor pressure and derivative.
   [~, des_sfc_dT] = icemodel.vapor.saturation_vapor_pressure( ...
      T_sfc, liqflag);

   % Bulk richardson stability function.
   stability = icemodel.surface.turbulence.bulk_richardson.stability_factor( ...
      T_sfc, tair, wspd, br_coefs);

   % Canonical non-conductive surface terms at the current linearization
   % state. The Robin solve handles conduction in the interior column system.
   [Qe, Qh] = ...
      icemodel.surface.turbulence.bulk_richardson.turbulent_heat_flux( ...
      T_sfc, tair, wspd, psfc, ea_atm, De, br_coefs, ro_air_Lv, liqflag, NaN);
   Qa = icemodel.surface.advective_heat_flux(ppt, tppt, cv_liq);
   Qsn = icemodel.surface.net_shortwave_radiation(Qsi, albedo, chi);
   [Qln, dQln_dT] = icemodel.surface.net_longwave_radiation(T_sfc, Qli);
   Q_sfc = icemodel.surface.evaluate_surface_energy_balance( ...
      Qsn, Qln, Qh, Qe, 0.0, Qa, 0.0);

   %%% Linearizations

   % Net longwave radiation.
   Fp_Qln = dQln_dT;

   % Sensible heat flux.
   Fp_Qh = -cv_air * De * stability;

   % Latent heat flux (linearization of es_sfc around T_sfc)
   % es(T) ≈ es + des_dT * (T - T_sfc) = (es - des_dT * T_sfc) + des_dT * T
   Fp_Qe = -ro_air_Lv * De * epsilon / psfc * stability * des_sfc_dT;

   % Total linearization.
   Fp = Fp_Qln + Fp_Qh + Fp_Qe;
   Fc = Q_sfc - Fp * T_sfc;
end
