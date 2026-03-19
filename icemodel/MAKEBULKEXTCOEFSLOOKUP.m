function lookup = MAKEBULKEXTCOEFSLOOKUP(dz_spect, spect_N, spect_S, solardwavl, varargin)
   %MAKEBULKEXTCOEFSLOOKUP Precompute bulk-extinction coefficients on a rho grid.
   %
   %  lookup = MAKEBULKEXTCOEFSLOOKUP(dz_spect, spect_N, spect_S, solardwavl)
   %  lookup = MAKEBULKEXTCOEFSLOOKUP(..., rho_grid)
   %
   % The lookup is approximate because it quantizes the snow density to the
   % supplied rho_grid before evaluating the spectral transform.
   %
   %#codegen

   % Default to an integer-density lookup over the physically relevant
   % range used by the current spectral model.
   if nargin > 4
      rho_grid = varargin{1};
   else
      rho_grid = transpose(300:917);
   end
   rho_grid = rho_grid(:);

   % Build one bulk-extinction profile per quantized density.
   n_bulk = size(spect_N, 1) + 2;
   table = zeros(n_bulk, numel(rho_grid));
   for i = 1:numel(rho_grid)
      table(:, i) = BULKEXTCOEFS(dz_spect, ...
         rho_grid(i) * ones(size(spect_N, 1), 1), ...
         spect_N, spect_S, solardwavl);
   end

   % Record the lookup metadata needed by the fast gather path.
   lookup = struct();
   lookup.rho_grid = rho_grid;
   lookup.bulkcoefs = table;
   lookup.rho_min = rho_grid(1);
   lookup.rho_max = rho_grid(end);
   lookup.rho_step = 0.0;
   lookup.is_uniform = false;
   if numel(rho_grid) > 1
      drho = diff(rho_grid);
      lookup.rho_step = drho(1);
      lookup.is_uniform = all(abs(drho - lookup.rho_step) <= 10 * eps(max(abs(rho_grid))));
   end
end
