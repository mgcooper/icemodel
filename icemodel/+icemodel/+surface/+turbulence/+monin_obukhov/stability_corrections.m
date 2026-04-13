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

   psi_m0_stable = icemodel.surface.turbulence.monin_obukhov.psi_holtslag(z0m ./ L);
   psi_mz_stable = icemodel.surface.turbulence.monin_obukhov.psi_holtslag(z_u ./ L);
   psi_h0_stable = icemodel.surface.turbulence.monin_obukhov.psi_holtslag(z0h ./ L);
   psi_hz_stable = icemodel.surface.turbulence.monin_obukhov.psi_holtslag(z_t ./ L);
   psi_q0_stable = icemodel.surface.turbulence.monin_obukhov.psi_holtslag(z0q ./ L);
   psi_qz_stable = icemodel.surface.turbulence.monin_obukhov.psi_holtslag(z_q ./ L);

   psi_m0_unstable = icemodel.surface.turbulence.monin_obukhov.psi_m_paulson(z0m ./ L);
   psi_mz_unstable = icemodel.surface.turbulence.monin_obukhov.psi_m_paulson(z_u ./ L);
   psi_h0_unstable = icemodel.surface.turbulence.monin_obukhov.psi_h_paulson(z0h ./ L);
   psi_hz_unstable = icemodel.surface.turbulence.monin_obukhov.psi_h_paulson(z_t ./ L);
   psi_q0_unstable = icemodel.surface.turbulence.monin_obukhov.psi_h_paulson(z0q ./ L);
   psi_qz_unstable = icemodel.surface.turbulence.monin_obukhov.psi_h_paulson(z_q ./ L);

   if isscalar(is_stable)
      if is_stable
         psi_m0 = psi_m0_stable;
         psi_mz = psi_mz_stable;
         psi_h0 = psi_h0_stable;
         psi_hz = psi_hz_stable;
         psi_q0 = psi_q0_stable;
         psi_qz = psi_qz_stable;
      else
         psi_m0 = psi_m0_unstable;
         psi_mz = psi_mz_unstable;
         psi_h0 = psi_h0_unstable;
         psi_hz = psi_hz_unstable;
         psi_q0 = psi_q0_unstable;
         psi_qz = psi_qz_unstable;
      end
      return
   end

   psi_m0 = psi_m0_unstable;
   psi_mz = psi_mz_unstable;
   psi_h0 = psi_h0_unstable;
   psi_hz = psi_hz_unstable;
   psi_q0 = psi_q0_unstable;
   psi_qz = psi_qz_unstable;

   psi_m0(is_stable) = psi_m0_stable(is_stable);
   psi_mz(is_stable) = psi_mz_stable(is_stable);
   psi_h0(is_stable) = psi_h0_stable(is_stable);
   psi_hz(is_stable) = psi_hz_stable(is_stable);
   psi_q0(is_stable) = psi_q0_stable(is_stable);
   psi_qz(is_stable) = psi_qz_stable(is_stable);
end
