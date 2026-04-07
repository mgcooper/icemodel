function [psi_m0, psi_mz, psi_h0, psi_hz, psi_q0, psi_qz] = ...
      stability_corrections(L, z0m, z_u, z0h, z_t, z0q, z_q, is_stable)
   %STABILITY_CORRECTIONS Return Monin-Obukhov profile corrections.
   %
   % Stable conditions use the Holtslag and de Bruin correction for both
   % momentum and scalars. Unstable conditions use Paulson/Dyer with separate
   % forms for momentum and scalars. The returned values are evaluated both at
   % the observation heights and at the roughness lengths because the flux
   % denominator uses psi(z/L) - psi(z0/L).
   %
   %#codegen

   if is_stable
      psi_m0 = icemodel.surface.turbulence.monin_obukhov.psi_holtslag(z0m / L);
      psi_mz = icemodel.surface.turbulence.monin_obukhov.psi_holtslag(z_u / L);
      psi_h0 = icemodel.surface.turbulence.monin_obukhov.psi_holtslag(z0h / L);
      psi_hz = icemodel.surface.turbulence.monin_obukhov.psi_holtslag(z_t / L);
      psi_q0 = icemodel.surface.turbulence.monin_obukhov.psi_holtslag(z0q / L);
      psi_qz = icemodel.surface.turbulence.monin_obukhov.psi_holtslag(z_q / L);
   else
      psi_m0 = icemodel.surface.turbulence.monin_obukhov.psi_m_paulson(z0m / L);
      psi_mz = icemodel.surface.turbulence.monin_obukhov.psi_m_paulson(z_u / L);
      psi_h0 = icemodel.surface.turbulence.monin_obukhov.psi_h_paulson(z0h / L);
      psi_hz = icemodel.surface.turbulence.monin_obukhov.psi_h_paulson(z_t / L);
      psi_q0 = icemodel.surface.turbulence.monin_obukhov.psi_h_paulson(z0q / L);
      psi_qz = icemodel.surface.turbulence.monin_obukhov.psi_h_paulson(z_q / L);
   end
end
