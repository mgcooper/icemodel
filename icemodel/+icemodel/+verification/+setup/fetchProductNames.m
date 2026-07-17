function names = fetchProductNames(registry)
   %FETCHPRODUCTNAMES Return ordered product selectors from a fetch registry.
   %
   %  names = icemodel.verification.setup.fetchProductNames(registry)
   %
   % Fetch registries store product selectors in column one and DOI provenance
   % in column two. This helper keeps that shared representation out of family
   % adapters while preserving the registry order used by their public defaults.

   arguments
      registry (:, 2) string
   end

   names = reshape(registry(:, 1), 1, []);
end
