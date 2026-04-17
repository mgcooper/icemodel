function [qext, g, coalbedo, wavel, dwavel] = get_scattering_coefficients( ...
      mie_table, opts)
   %get_scattering_coefficients Extract spectral scattering coefficients from
   %a pre-loaded Mie scattering matrix and the wavelength grid.
   %
   % [qext, g, coalbedo, wavel, dwavel] = ...
   %    icemodel.radiation.get_scattering_coefficients(mie_table, opts)
   %
   % MIE_TABLE is the raw Mie scattering matrix already loaded from disk (e.g.
   % by icemodel.radiation.load_spectral_tables). The matrix stores one row
   % block per optical property:
   %   rows 1:nradii                 -> g
   %   rows nradii+1:2*nradii        -> qext
   %   rows 2*nradii+1:3*nradii      -> coalbedo
   %   rows 3*nradii+1:4*nradii      -> wavel
   %
   % The wavelength grid is repeated for every grain radius in the table, so
   % only the canonical wavelength vector is returned here.
   %
   %#codegen

   % Extract the single-scattering property tables for the configured number of
   % radii and wavelengths.
   g = mie_table(1:opts.nradii, 1:opts.nwavel);
   qext = mie_table(opts.nradii + 1:2 * opts.nradii, 1:opts.nwavel);
   coalbedo = mie_table(2 * opts.nradii + 1:3 * opts.nradii, 1:opts.nwavel);
   wavel = mie_table(3 * opts.nradii + 1, 1:opts.nwavel);

   % Build the quadrature weights from the fixed 118-band wavelength grid. The
   % end bands use one-sided widths, doubled so the discrete quadrature still
   % spans the full first and last wavelength intervals.
   dwavel = zeros(1, opts.nwavel);
   dwavel(1) = 2.0 * (wavel(2) - wavel(1));
   dwavel(2:end - 1) = (wavel(3:end) - wavel(1:end - 2)) / 2.0;
   dwavel(end) = 2.0 * (wavel(end) - wavel(end - 1));
end
