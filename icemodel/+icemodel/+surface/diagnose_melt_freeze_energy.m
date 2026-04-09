function [Qm, Qf] = diagnose_melt_freeze_energy(T_sfc, chi, albedo, Qsi, ...
      Qli, Qle, Qh, Qe, Qc, Qa, tair, wspd, psfc, ppt, tppt, De, ea_atm, ...
      br_coefs, roL, T, k_eff, dz, ro_sfc, snow_depth, opts)
   %DIAGNOSE_MELT_FREEZE_ENERGY Diagnose surplus/deficit energy relative to Tf.
   %
   %  [Qm, Qf] = icemodel.surface.diagnose_melt_freeze_energy(...)
   %
   %  Qm > 0 is the energy surplus available for melting ice.
   %  Qf > 0 is the energy deficit required to warm the surface to Tf.
   %
   % The flux terms (Qle, Qh, Qe, Qc, Qa) passed in must be evaluated at
   % T_sfc_phys = min(T_sfc, Tf), i.e., after the physical-surface-temperature
   % cap is applied.
   %
   % When T_sfc < Tf, the freezing deficit Qf is diagnosed by re-evaluating
   % the full SEB at Tf rather than at the current T_sfc.
   %
   % See also: icemodel.surface.surface_energy_balance_residual,
   %           icemodel.surface.evaluate_surface_energy_balance
   %
   %#codegen

   persistent Tf emiss
   if isempty(Tf)
      Tf = icemodel.physicalConstant('Tf');
      emiss = icemodel.parameterLookup('emiss');
   end

   Qm = 0.0;
   Qf = 0.0;

   if T_sfc >= Tf
      % Compute melt energy: the balance at the melt cap with no melt term.
      Qsn = chi * Qsi * (1.0 - albedo);
      Qln = emiss * Qli + Qle;
      Qm = icemodel.surface.evaluate_surface_energy_balance( ...
         Qsn, Qln, Qh, Qe, Qc, Qa, 0.0);
   else
      % Compute energy needed to reach melt temp (energy deficit) by
      % re-evaluating the full SEB at Tf. Qc must be recomputed at Tf here.
      Qf = -(icemodel.surface.surface_energy_balance_residual(Tf, tair, Qsi, ...
         Qli, albedo, wspd, ppt, tppt, psfc, De, ea_atm, br_coefs, roL, true, chi, ...
         icemodel.surface.conductive_heat_flux(k_eff, T, dz, Tf), ...
         ro_sfc, snow_depth, opts));
   end
end
