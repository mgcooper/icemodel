function [radii, mie, solar_spectrum, kabs, kice] = load_spectral_tables(opts)
   %LOAD_SPECTRAL_TABLES Load raw spectral input tables from disk.
   %
   %  [radii, mie, solar_spectrum, kabs, kice] = ...
   %     icemodel.radiation.load_spectral_tables(opts)
   %
   % Returns the unprocessed spectral lookup tables used by the two-stream
   % radiation model. Useful for ad hoc inspection or re-processing outside
   % the standard initialization path.
   %
   % For production use, call icemodel.radiation.initialize_spectral_model,
   % which loads these tables internally and returns the processed optical
   % coefficients and spectral mesh ready for the solver.
   %
   %#codegen

   % The grain radii available to the two-stream spectral model [mm].
   radii = [ ...
      0.040, 0.050, 0.065, 0.080, 0.100, ...
      0.120, 0.140, 0.170, 0.200, 0.240, ...
      0.290, 0.350, 0.420, 0.500, 0.570, ...
      0.660, 0.760, 0.870, 1.000, 1.100, ...
      1.250, 1.400, 1.600, 1.800, 2.000, ...
      2.250, 2.500, 2.750, 3.000, 3.500, ...
      4.000, 4.500, 5.000, 5.500, 6.000 ];

   % Load the pre-defined Mie scattering values.
   mie = load(fullfile(opts.pathinput, 'spectral', 'mie.mat'), 'mie');
   mie = mie.mie;

   % Load the prototype spectral irradiance profile.
   solar_spectrum = load( ...
      fullfile(opts.pathinput, 'spectral', 'solar.mat'), 'solar');
   solar_spectrum = solar_spectrum.solar;

   % Load the optional user-defined absorption coefficients.
   if opts.kabs_user
      kabs = load(fullfile(opts.pathinput, 'spectral', 'kabs.mat'), 'kabs');
      kice = load(fullfile(opts.pathinput, 'spectral', 'kice.mat'), 'kice');
      kabs = kabs.kabs;
      kice = kice.kice;
   else
      kabs = [];
      kice = [];
   end
end
