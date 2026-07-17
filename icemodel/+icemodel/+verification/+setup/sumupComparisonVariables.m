function variables = sumupComparisonVariables(observations)
   %SUMUPCOMPARISONVARIABLES List nonempty SUMup observation groups.
   %
   %  variables = ...
   %     icemodel.verification.setup.sumupComparisonVariables(observations)

   % Keep the canonical comparison order while omitting groups that were not
   % staged for the selected point and observation window.
   candidates = ["density"; "subsurface_temperature"; "smb"];
   present = false(numel(candidates), 1);
   for k = 1:numel(candidates)
      present(k) = isfield(observations, candidates(k)) ...
         && ~isempty(observations.(candidates(k)));
   end
   variables = candidates(present);
end
