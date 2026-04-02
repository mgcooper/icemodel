function [Qm, Qf] = MFENERGY(albedo, Qsi, Qli, Qle, Qh, Qe, Qc, Qa, Ts, ...
      Tf, Ta, wspd, ppt, tppt, De, ea, roL, Pa, k_eff, T, dz, scoef, ...
      chi, ro_sfc, snow_depth, opts)
   %MFENERGY Compute the energy available for melting or freezing.
   %
   %  Qm > 0 is the energy surplus available for melting ice.
   %  Qf > 0 is the energy deficit required to warm the surface to Tf.
   %
   % The fluxes passed into this function must be evaluated using
   % Ts = min(Ts, Tf). When Ts < Tf, the freezing deficit is diagnosed by
   % re-evaluating the canonical surface energy balance at Ts = Tf.
   %
   %#codegen

   persistent emiss
   if isempty(emiss)
      emiss = icemodel.parameterLookup('emiss');
   end

   Qm = 0.0;
   Qf = 0.0;

   if Ts >= Tf
      % Compute melt energy
      Qm = ENBAL(albedo, emiss, chi, Qsi, Qli, Qle, Qh, Qe, Qc, Qa, 0.0);
   else
      % Compute energy needed to reach melt temp (energy deficit)
      [~, ~, ~, ~, ~, ~, balance] = ...
         icemodel.surface.surface_energy_balance_terms(Tf, Ta, Qsi, Qli, ...
         albedo, wspd, ppt, tppt, Pa, De, ea, T, k_eff, dz, roL, chi, ...
         scoef, true, ro_sfc, snow_depth, opts);
      Qf = -balance;
   end
end
