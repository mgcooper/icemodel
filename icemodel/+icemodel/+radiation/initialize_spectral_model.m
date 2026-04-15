function [I0, dz, z_nodes, z_edges, tau_N, tau_S, solar_dwavel, ...
      k_bulk_lookup, r_eff, qext, g, coalbedo, kabs, kice, wavel, radii] = ...
      initialize_spectral_model(opts)
   %initialize_spectral_model Initialize spectral geometry and coefficients.
   %
   % [I0, dz, z_nodes, z_edges, tau_N, tau_S, solar_dwavel, ...
   %    k_bulk_lookup, r_eff] ...
   %    = icemodel.radiation.initialize_spectral_model(opts)
   %
   % [I0, dz, z_nodes, z_edges, tau_N, tau_S, solar_dwavel, ...
   %    k_bulk_lookup, r_eff, qext, g, coalbedo, kabs, kice, wavel, radii] ...
   %     = icemodel.radiation.initialize_spectral_model(opts)
   %
   % z_nodes are the centers of the spectral control volumes. z_edges includes
   % the top surface of the top grid cell and the bottom surface of the bottom
   % grid cell. tau_N and tau_S are the precomputed edge-based optical-depth
   % terms used by the exact bulk-extinction coefficient transform.
   %
   % Spectral input tables (Mie scattering matrix, solar spectrum, and optional
   % absorption profiles) are loaded via icemodel.radiation.load_spectral_tables.
   %
   % The additional optical-property outputs are returned so a future grain-size
   % evolution model can refresh k_ext through
   % `icemodel.radiation.update_extinction_coefficients` without reloading the
   % Mie tables every timestep. r_eff is the configured optical grain radius
   % from the Mie lookup table, which `icemodel.column.initialize_column_state`
   % expands onto the thermal mesh.
   %
   %#codegen

   % Build the spectral mesh.
   [dz, ~, z_nodes, z_edges] = ...
      icemodel.column.control_volume_mesh(opts.z0_spectral, opts.dz_spectral);

   % Load raw spectral input tables.
   [radii, mie_table, solar_spectrum, kabs, kice] = ...
      icemodel.radiation.load_spectral_tables(opts);

   % Load the single-scattering property tables and the fixed wavelength grid.
   [qext, g, coalbedo, wavel, dwavel] = ...
      icemodel.radiation.get_scattering_coefficients(mie_table, opts);

   % Integrate the spectral irradiance profile over the model wavelength grid.
   [I0, ~, solar_dwavel] = ...
      icemodel.radiation.get_solar_spectrum(solar_spectrum, wavel, dwavel);

   % Build k_ext and its edge-based optical-depth terms for the configured
   % optical grain radius, and the bulk-extinction lookup table if requested.
   [tau_N, tau_S, k_bulk_lookup] = ...
      icemodel.radiation.update_extinction_coefficients( ...
      qext, g, coalbedo, kabs, kice, wavel, radii, opts.i_grainradius, ...
      z_edges, dz, solar_dwavel, opts.lookup_k_bulk);

   % Return the initial optical grain radius [mm] from the lookup table.
   % Note: this is the optically equivalent radius from the Mie tables. It's
   % currently used to initialize the thermal grain radius tracked by
   % icemodel.column.vapor_mass_transfer, but these concepts are not
   % identical; coupling them is future
   % work (see icemodel.radiation.update_extinction_coefficients).
   r_eff = radii(opts.i_grainradius);
end
