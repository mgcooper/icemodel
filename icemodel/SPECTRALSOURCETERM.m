function [Sc, chi] = SPECTRALSOURCETERM(Qsi, albedo, Q0, Q, dz_spect, ...
      dz_therm, template)
   %SPECTRALSOURCETERM Convert net spectral flux to absorbed source term.
   %
   %  [Sc, chi] = SPECTRALSOURCETERM(Qsi, albedo, Q0, Q, dz_spect, ...
   %     dz_therm, template)
   %
   % This helper reproduces the source-term construction embedded in
   % UPDATEEXTCOEFS so alternative spectral implementations can share the
   % same absorbed-flux and SEB-partition logic.
   %
   %#codegen

   % Collapse the spectral-grid net-flux differences back onto the thermal
   % grid used by the main column solve.
   dQ = Q(1:end-1) - Q(2:end);
   dQ = transpose(sum(reshape(dQ, dz_therm(1) / dz_spect(1), []), 1));
   dQ = vertcat(dQ, zeros((sum(dz_therm) - sum(dz_spect)) / dz_therm(1), 1));

   % Preserve the established top-layer chi rule used by the SEB coupling.
   if albedo > 0.65
      chi = 0.9;
   else
      chi = dQ(1) / sum(dQ);
   end

   % Convert absorbed flux to a volumetric source term and reshape it to the
   % caller's template when the spectral helper is fed row-vector state.
   Sc = -(1.0 - chi) * Qsi / Q0 * dQ ./ dz_therm;
   if ~isequal(size(Sc), size(template))
      Sc = reshape(Sc, size(template));
   end
end
