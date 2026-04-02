function [Qe, Qh, Qc, Qm, Qf, Qbal] = SEBFLUX(T, xTs, Ta, Qsi, Qli, ...
      albedo, wspd, ppt, tppt, Pa, De, ea, Tf, k_eff, dz, roL, chi, ...
      scoef, liqflag, ro_sfc, snow_depth, opts)
   %SEBFLUX using the new top node temperature, compute a new surface flux
   %
   %
   % ea  - atmospheric vapor pressure from relative humidity data.
   % S   - stability function.
   % es  - surface saturation water vapor pressure.
   % Qe  - turbulent latent heat flux.
   % Qh  - turbulent sensible heat flux.
   % Qle - surface emitted longwave heat flux.
   % Qc  - conductive heat flux into the surface.
   % Qm  - the energy flux available for melting or freezing.
   % bal - the energy balance residual.
   %
   %#codegen

   persistent emiss
   if isempty(emiss)
      emiss = icemodel.parameterLookup('emiss');
   end

   Ts = MELTTEMP(xTs,Tf);

   [Ts, Qe, Qh, Qc, Qa, Qle, balance] = ...
      icemodel.surface.surface_energy_balance_terms(Ts, Ta, Qsi, ...
      Qli, albedo, wspd, ppt, tppt, Pa, De, ea, T, k_eff, dz, ...
      roL, chi, scoef, liqflag, ro_sfc, snow_depth, opts);

   [Qm,Qf] = MFENERGY(albedo,Qsi,Qli,Qle,Qh,Qe,Qc,Qa,xTs,Tf, ... % Note: xTs
      Ta,wspd,ppt,tppt,De,ea,roL,Pa,k_eff,T,dz,scoef,chi,ro_sfc, ...
      snow_depth,opts);

   Qbal = ENBAL(albedo,emiss,chi,Qsi,Qli,Qle,Qh,Qe,Qc,Qa,Qm);
end
