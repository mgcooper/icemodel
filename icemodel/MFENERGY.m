function [Qm, Qf] = MFENERGY(T_sfc, chi, albedo, Qsi, Qli, Qle, Qh, Qe, Qc, ...
      Qa, tair, wspd, psfc, ppt, tppt, ea_atm, br_coefs, De, T, k_eff, dz, ...
      ro_sfc, snow_depth, roL, opts)
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

   persistent Tf
   if isempty(Tf)
      Tf = icemodel.physicalConstant('Tf');
   end

   Qm = 0.0;
   Qf = 0.0;

   if T_sfc >= Tf
      % Compute melt energy
      Qm = ENBAL(chi, albedo, Qsi, Qli, Qle, Qh, Qe, Qc, Qa, 0.0);
   else
      % Compute energy needed to reach melt temp (energy deficit)
      [~, ~, ~, ~, ~, ~, balance] = ...
         icemodel.surface.surface_energy_balance_terms(Tf, tair, Qsi, Qli, ...
         albedo, wspd, ppt, tppt, psfc, De, ea_atm, T, k_eff, dz, roL, chi, ...
         br_coefs, true, ro_sfc, snow_depth, opts);
      Qf = -balance;
   end
end
