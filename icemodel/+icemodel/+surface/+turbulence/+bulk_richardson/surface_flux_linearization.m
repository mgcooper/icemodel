function [Fc, Fp] = surface_flux_linearization(tair, Qsi, Qli, albedo, wspd, ...
      ppt, tppt, psfc, De, ea_atm, roL, br_coefs, chi, T_sfc, liqflag)
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
   % See also: icemodel.surface.turbulence.bulk_richardson.surface_fluxes
   %
   %#codegen

   persistent cv_air cv_liq emiss SB epsilon
   if isempty(cv_air)
      [cv_air, cv_liq, SB, epsilon] = icemodel.physicalConstant( ...
         'cv_air', 'cv_liq', 'SB', 'epsilon');
      emiss = icemodel.parameterLookup('emiss');
   end

   % Surface saturation vapor pressure and derivative from icemodel.vapor.vappress
   [es_sfc, des_sfc_dT] = icemodel.vapor.vappress(T_sfc, liqflag);

   % Bulk richardson stability function.
   stability = icemodel.surface.turbulence.bulk_richardson.stability_factor( ...
      T_sfc, tair, wspd, br_coefs);

   %%% Linearizations

   % Emitted longwave (note: -Qle = -emiss * SB * T_sfc ^ 4)
   Fc_Qle = 3 * emiss * SB * T_sfc ^ 4;
   Fp_Qle = -4 / 3 * Fc_Qle / T_sfc;

   % Sensible heat flux
   Fc_Qh = cv_air * De * stability * tair;
   Fp_Qh = -cv_air * De * stability;

   % Latent heat flux (linearization of es_sfc around T_sfc)
   % es(T) ≈ es + des_dT * (T - T_sfc) = (es - des_dT * T_sfc) + des_dT * T
   Fc_Qe = roL * De * epsilon / psfc * stability ...
      * (ea_atm - es_sfc + des_sfc_dT * T_sfc);
   Fp_Qe = -roL * De * epsilon / psfc * stability * des_sfc_dT;

   % Precipitation-advected heat:
   Fc_Qa = QADVECT(ppt, tppt, cv_liq);
   Fp_Qa = 0;

   % Combine net sw, incoming lw, and precipitation-advected heat:
   Fc = Fc_Qle + Fc_Qh + Fc_Qe + emiss * Qli + chi * Qsi * (1.0 - albedo) ...
      + Fc_Qa;
   Fp = Fp_Qle + Fp_Qh + Fp_Qe + Fp_Qa;
end
