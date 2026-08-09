function assertNotEvaluationDestination(destinations, evaluation_roots)
   %ASSERTNOTEVALUATIONDESTINATION Refuse reconstruction writes under eval data.
   %
   %  icemodel.forcing.reconstruct.assertNotEvaluationDestination( ...
   %     destinations, evaluation_roots)
   %
   % Reconstruction may read evaluation observations, but its persistence
   % boundaries must never create, replace, or remove files in an evaluation
   % tree. Canonical containment closes relative-path and symlink aliases.

   destinations = reshape(string(destinations), [], 1);
   evaluation_roots = unique(reshape(string(evaluation_roots), [], 1));
   evaluation_roots(ismissing(evaluation_roots) ...
      | strlength(evaluation_roots) == 0) = [];

   % Test every destination against every protected root before a caller
   % creates directories or opens a file.
   for destination = destinations'
      for root = evaluation_roots'
         if icemodel.isPathInside(destination, root)
            error(['icemodel:reconstruct:' ...
               'assertNotEvaluationDestination:protectedPath'], ...
               'Reconstruction cannot write beneath evaluation root %s: %s', ...
               root, destination)
         end
      end
   end
end
