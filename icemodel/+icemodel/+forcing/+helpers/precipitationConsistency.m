function valid = precipitationConsistency(ppt, rainf, snowf)
   %PRECIPITATIONCONSISTENCY Check finite, nonnegative phase mass balance.
   %
   %  valid = icemodel.forcing.helpers.precipitationConsistency( ...
   %     ppt, rainf, snowf)
   %
   % Returns one logical value per sample. The tolerance is the single source
   % used by reconstruction and verification artifact audits.

   arguments
      ppt (:, 1) double
      rainf (:, 1) double
      snowf (:, 1) double
   end

   % Refuse axes that cannot represent a samplewise partition.
   if numel(rainf) ~= numel(ppt) || numel(snowf) ~= numel(ppt)
      error('icemodel:forcing:precipitationConsistency:sizeMismatch', ...
         'ppt, rainf, and snowf must share one sample axis');
   end

   % Admit only finite, nonnegative components whose sum preserves total mass.
   tolerance = max(1e-12, 1e-9 .* abs(ppt));
   valid = isfinite(ppt) & isfinite(rainf) & isfinite(snowf) ...
      & ppt >= 0 & rainf >= 0 & snowf >= 0 ...
      & abs(rainf + snowf - ppt) <= tolerance;
end
