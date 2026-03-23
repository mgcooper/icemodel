function k_ext = SPECTEXTCOEFS(qext, g, coalbedo, radii, iradius)
   %SPECTEXTCOEFS Compute spectral extinction coefficients for one grain radius.
   %
   %  k_ext = SPECTEXTCOEFS(qext, g, coalbedo, radii, iradius)
   %
   % The optical-property tables are loaded once and indexed here by IRADIUS.
   % Integer IRADIUS values select one tabulated grain size directly. A
   % fractional IRADIUS linearly interpolates between the two neighboring table
   % rows, which is the scaffold needed if a future grain-size model maps its
   % evolving optical grain size onto the Mie table.
   %
   %#codegen

   n_radii = numel(radii);
   iradius = min(max(iradius, 1.0), n_radii);
   i0 = floor(iradius);
   i1 = ceil(iradius);
   w1 = iradius - i0;
   w0 = 1.0 - w1;

   % Convert grain radius from mm to um.
   r_snow = (w0 * radii(i0) + w1 * radii(i1)) / 1000.0;

   % Compute the optical properties for this grain radius.
   g_r = w0 * g(i0, :) + w1 * g(i1, :);
   qext_r = w0 * qext(i0, :) + w1 * qext(i1, :);
   coalbedo_r = w0 * coalbedo(i0, :) + w1 * coalbedo(i1, :);

   % Compute the spectral extinction coefficients.
   sigma_e = (3.0 / 4.0) * qext_r / r_snow;
   k_ext = sigma_e .* sqrt(coalbedo_r - coalbedo_r .* g_r ...
      + coalbedo_r .^ 2 .* g_r);
end
