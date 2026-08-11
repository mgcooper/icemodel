function assertNotEvaluationDestination(destinations, evaluation_roots)
   %ASSERTNOTEVALUATIONDESTINATION Refuse reconstruction writes under eval data.
   %
   %  icemodel.forcing.reconstruct.assertNotEvaluationDestination( ...
   %     destinations, evaluation_roots)
   %
   % Reconstruction may read evaluation observations. It must never create,
   % replace, or remove a file in an evaluation tree. The check compares
   % canonical paths, so a relative path or a symlink cannot get past it.

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
