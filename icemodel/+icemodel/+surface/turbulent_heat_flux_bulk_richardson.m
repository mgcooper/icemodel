function [Qe, Qh, diag] = turbulent_heat_flux_bulk_richardson(Ta, Ts, ...
      wspd, Pa, De, ea_atm, cv_air, roL, scoef, liqflag, z0_bulk)
   %TURBULENT_HEAT_FLUX_BULK_RICHARDSON Evaluate the explicit THF scheme.
   %
   %  [Qe, Qh] = icemodel.surface.turbulent_heat_flux_bulk_richardson(...)
   %  [Qe, Qh, diag] = icemodel.surface.turbulent_heat_flux_bulk_richardson(...)
   %
   % This is the existing bulk-Richardson surface-transfer scheme used by the
   % production SEB path. The formulation uses one aerodynamic roughness, one
   % exchange coefficient, and the Louis/Liston stability factor.
   %
   %#codegen

   persistent epsilon
   if isempty(epsilon)
      epsilon = icemodel.physicalConstant('epsilon');
   end

   es_sfc = VAPPRESS(Ts, liqflag);
   St = STABLEFN(Ta, Ts, wspd, scoef);
   Qe = LATENT(De, St, ea_atm, es_sfc, roL, epsilon, Pa);
   Qh = SENSIBLE(De, St, Ta, Ts, cv_air);

   if nargout > 2
      diag = bulk_richardson_diag(St, es_sfc, z0_bulk);
   end
end

function diag = bulk_richardson_diag(stability_factor, es_sfc, z0_bulk)
   %BULK_RICHARDSON_DIAG Assemble optional diagnostics for test/debug use.

   diag = struct( ...
      'scheme', 'bulk_richardson', ...
      'stability_factor', stability_factor, ...
      'es_sfc', es_sfc, ...
      'z0m', z0_bulk, ...
      'z0h', NaN, ...
      'z0q', NaN, ...
      'u_star', NaN, ...
      'L', NaN, ...
      'Re', NaN);
end
