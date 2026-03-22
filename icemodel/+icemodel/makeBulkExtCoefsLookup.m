function lookup = makeBulkExtCoefsLookup(dz_spect, tau_N, ...
      tau_S, solar_dwavel, varargin)
   %MAKEBULKEXTCOEFSLOOKUP Precompute bulk-extinction coefficients on a density grid.
   %
   %  lookup = icemodel.makeBulkExtCoefsLookup(dz_spect, ...
   %     tau_N, tau_S, solar_dwavel)
   %  lookup = icemodel.makeBulkExtCoefsLookup(..., ro_grid)
   %
   % The lookup is approximate because it quantizes snow density to the
   % supplied RO_GRID before evaluating the spectral transform.

   %#codegen

   % Default to an integer-density lookup over the current supported range.
   if nargin > 4
      ro_grid = varargin{1};
   else
      ro_grid = transpose(300:917);
   end
   ro_grid = ro_grid(:);

   % Build one bulk-extinction profile per quantized density.
   n_bulk = size(tau_N, 1) + 2;
   table = zeros(n_bulk, numel(ro_grid));
   for i = 1:numel(ro_grid)
      table(:, i) = BULKEXTCOEFS(dz_spect, ...
         ro_grid(i) * ones(size(tau_N, 1), 1), ...
         tau_N, tau_S, solar_dwavel);
   end

   % Record the lookup metadata needed by the fast gather path.
   lookup = struct();
   lookup.ro_grid = ro_grid;
   lookup.bulkcoefs = table;
   lookup.ro_min = ro_grid(1);
   lookup.ro_max = ro_grid(end);
   lookup.ro_step = 0.0;
   lookup.is_uniform = false;
   if numel(ro_grid) > 1
      dro = diff(ro_grid);
      lookup.ro_step = dro(1);
      lookup.is_uniform = ...
         all(abs(dro - lookup.ro_step) <= 10 * eps(max(abs(ro_grid))));
   end
end
