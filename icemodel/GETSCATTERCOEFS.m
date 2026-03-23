function [qext, g, coalbedo, wavel, dwavel, radii] = GETSCATTERCOEFS(opts)
   %GETSCATTERCOEFS Load the spectral scattering tables and wavelength grid.
   %
   %  [qext, g, coalbedo, wavel, dwavel, radii] = GETSCATTERCOEFS(opts)
   %
   % The Mie table stores one row block per optical property:
   %   rows 1:nradii                 -> g
   %   rows nradii+1:2*nradii        -> qext
   %   rows 2*nradii+1:3*nradii      -> coalbedo
   %   rows 3*nradii+1:4*nradii      -> wavel
   %
   % The wavelength grid is repeated for every grain radius in the table, so
   % only the canonical wavelength vector is returned here.
   %
   %#codegen

   % The grain radii available to the two-stream spectral model.
   radii = [ ...
      0.040, 0.050, 0.065, 0.080, 0.100, ...
      0.120, 0.140, 0.170, 0.200, 0.240, ...
      0.290, 0.350, 0.420, 0.500, 0.570, ...
      0.660, 0.760, 0.870, 1.000, 1.100, ...
      1.250, 1.400, 1.600, 1.800, 2.000, ...
      2.250, 2.500, 2.750, 3.000, 3.500, ...
      4.000, 4.500, 5.000, 5.500, 6.000 ];

   % Load the pre-defined Mie scattering values once.
   scattercoefs = load(fullfile(opts.pathinput, 'spectral', 'mie.mat'));
   scattercoefs = scattercoefs.mie;

   % Extract the single-scattering property tables for the configured number of
   % radii and wavelengths.
   g = scattercoefs(1:opts.nradii, 1:opts.nwavel);
   qext = scattercoefs(opts.nradii + 1:2 * opts.nradii, 1:opts.nwavel);
   coalbedo = scattercoefs(2 * opts.nradii + 1:3 * opts.nradii, 1:opts.nwavel);
   wavel = scattercoefs(3 * opts.nradii + 1, 1:opts.nwavel);

   % Build the quadrature weights from the fixed 118-band wavelength grid. The
   % end bands use one-sided widths, doubled so the discrete quadrature still
   % spans the full first and last wavelength intervals.
   dwavel = zeros(1, opts.nwavel);
   dwavel(1) = 2.0 * (wavel(2) - wavel(1));
   dwavel(2:end - 1) = (wavel(3:end) - wavel(1:end - 2)) / 2.0;
   dwavel(end) = 2.0 * (wavel(end) - wavel(end - 1));
end
