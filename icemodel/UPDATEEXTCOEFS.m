function [k_ext, tau_N, tau_S] = UPDATEEXTCOEFS(qext, g, coalbedo, radii, ...
      iradius, z_edges, ro_ice, wavel, kice, kabs)
   %UPDATEEXTCOEFS Update spectral extinction coefficients for one grain radius.
   %
   %  [k_ext, tau_N, tau_S] = UPDATEEXTCOEFS(qext, g, coalbedo, radii, ...
   %     iradius, z_edges, ro_ice, wavel, kice, kabs)
   %
   % This helper computes k_ext for the requested optical grain-radius index,
   % or interpolated index, applies the optional impurity scaling, then
   % precomputes tau_N and tau_S for the exact bulk-extinction transform. It is
   % the natural entry point if a future grain-growth model needs to refresh
   % k_ext during the timestep loop without reloading the optical tables.
   %
   %#codegen

   % Compute the spectral extinction coefficients.
   k_ext = SPECTEXTCOEFS(qext, g, coalbedo, radii, iradius);
   if xor(isempty(kice), isempty(kabs))
      error('kice and kabs must both be empty or both be provided')
   end
   if ~isempty(kice)
      k_ext = SCALESPECTEXTCOEFS(k_ext, wavel, kice, kabs);
   end

   % Compute specific optical depth (precomputed arguments to the exponential
   % function in BULKEXTCOEFS)
   tau_edges = -(z_edges * k_ext) / ro_ice;
   tau_N = tau_edges(1:end - 1, :);
   tau_S = tau_edges(2:end, :);
end
