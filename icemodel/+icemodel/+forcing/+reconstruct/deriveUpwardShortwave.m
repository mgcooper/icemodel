function [filled, provenance, audit] = deriveUpwardShortwave( ...
      filled, provenance, audit, codes)
   %DERIVEUPWARDSHORTWAVE Fill missing swu from final albedo and swd.
   %
   %  [filled, provenance, audit] = ...
   %     icemodel.forcing.reconstruct.deriveUpwardShortwave( ...
   %     filled, provenance, audit, codes)
   %
   % Role
   %  Enforce the policy ordering rule swu = albedo * swd only after both
   %  operand channels have completed reconstruction. Native finite swu is
   %  preserved; unresolved operands leave swu missing.
   %
   % See also: icemodel.forcing.reconstruct.fillPromiceStation,
   %  icemodel.forcing.reconstruct.provenanceCodes

   arguments
      filled timetable
      provenance timetable
      audit table
      codes (1, 1) struct
   end

   required = ["swu", "albedo", "swd"];
   if ~all(ismember(required, string(filled.Properties.VariableNames))) ...
         || ~ismember("swu", string(provenance.Properties.VariableNames))
      error('icemodel:reconstruct:deriveUpwardShortwave:missingChannel', ...
         'swu derivation requires swu, albedo, swd, and swu provenance');
   end

   % The algebraic relation is valid only where both final operands exist.
   derive = ~isfinite(filled.swu) & isfinite(filled.albedo) ...
      & isfinite(filled.swd);
   if ~any(derive)
      return
   end
   candidate = filled.albedo(derive) .* filled.swd(derive);
   valid = icemodel.forcing.reconstruct.physicalValidity( ...
      "swu", candidate, filled.Properties.RowTimes(derive), ...
      swd=filled.swd(derive));
   target = find(derive);
   target = target(valid);
   if isempty(target)
      return
   end

   filled.swu(target) = filled.albedo(target) .* filled.swd(target);
   provenance.swu(target) = codes.derived_shortwave;
   derived = false(height(filled), 1);
   derived(target) = true;
   rows = icemodel.forcing.reconstruct.auditSegments( ...
      filled.Properties.RowTimes, derived, "swu", ...
      "derived_shortwave", ...
      "swu = albedo * swd after final operand reconstruction");
   audit = [audit; cell2table(vertcat(rows{:}), ...
      'VariableNames', audit.Properties.VariableNames)];
end
