function [I0, dz, z_nodes, z_edges, tau_N, tau_S, solar_dwavel, ...
      k_bulk_lookup, qext, g, coalbedo, radii, wavel, dwavel, k_ext, ...
      kice, kabs] = EXTCOEFSINIT(opts, ro_ice)
   %EXTCOEFSINIT Initialize the spectral geometry and optical coefficients.
   %
   %  [I0, dz, z_nodes, z_edges, tau_N, tau_S, solar_dwavel, ...
   %     k_bulk_lookup] = EXTCOEFSINIT(opts, ro_ice)
   %  [I0, dz, z_nodes, z_edges, tau_N, tau_S, solar_dwavel, ...
   %     k_bulk_lookup, qext, g, coalbedo, radii, wavel, dwavel, ...
   %     k_ext, kice, kabs] = EXTCOEFSINIT(opts, ro_ice)
   %
   % z_nodes are the centers of the spectral control volumes. z_edges includes
   % the top surface of the top grid cell and the bottom surface of the bottom
   % grid cell. tau_N and tau_S are the precomputed edge-based optical-depth
   % terms used by the exact bulk-extinction transform.
   %
   % The additional optical-property outputs are returned so a future grain-size
   % evolution model can refresh k_ext through UPDATEEXTCOEFS without reloading
   % the Mie tables every timestep.
   %
   %#codegen

   % Build the spectral mesh.
   [dz, ~, z_nodes, z_edges] = CVMESH(opts.z0_spectral, opts.dz_spectral);

   % Load the single-scattering property tables and the fixed wavelength grid.
   [qext, g, coalbedo, wavel, dwavel, radii] = GETSCATTERCOEFS(opts);

   % Load the prototype spectral irradiance profile and integrate it on the
   % model wavelength grid.
   solar_spectrum = load(fullfile(opts.pathinput, 'spectral', 'solar.mat'));
   solar_spectrum = solar_spectrum.solar;
   [I0, ~, solar_dwavel] = GETSOLAR(solar_spectrum, wavel, dwavel);

   % Load the optional user-defined absorption profile.
   if opts.kabs_user
      kabs = load(fullfile(opts.pathinput, 'spectral', 'kabs.mat'));
      kice = load(fullfile(opts.pathinput, 'spectral', 'kice.mat'));
      kabs = kabs.kabs;
      kice = kice.kice;
   else
      kabs = [];
      kice = [];
   end

   % Build k_ext and its edge-based optical-depth terms for the configured
   % optical grain radius.
   [k_ext, tau_N, tau_S] = UPDATEEXTCOEFS(qext, g, coalbedo, radii, ...
      opts.i_grainradius, z_edges, ro_ice, wavel, kice, kabs);

   % Build the bulk-extinction lookup table when requested.
   if isfield(opts, 'lookup_k_bulk') && opts.lookup_k_bulk
      k_bulk_lookup = icemodel.makeBulkExtCoefsLookup( ...
         dz, tau_N, tau_S, solar_dwavel);
   else
      k_bulk_lookup = struct([]);
   end
end
