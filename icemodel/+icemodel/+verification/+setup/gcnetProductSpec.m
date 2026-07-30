function spec = gcnetProductSpec(product)
   %GCNETPRODUCTSPEC Return Vandecrux/GC-Net DOI metadata and file rules.
   %
   %  spec = icemodel.verification.setup.gcnetProductSpec("surface")
   %
   %  This is the shared metadata authority for GC-Net/Vandecrux fetch validation
   %  and inventory discovery. gcnetProductNames owns accepted selector order;
   %  keep DOI metadata and station suffixes here so discovery cannot drift.

   arguments
      product (1, 1) string ...
         {icemodel.verification.validators.mustBeGcnetProductSelection(product)}
   end

   % Each product maps to one Arctic Data DOI package, one EML metadata file,
   % and the station-specific NetCDF suffixes required by RetMIP staging.
   switch product
      case "surface"
         spec = struct( ...
            'product', product, ...
            'doi', "10.18739/A2HM52K87", ...
            'metadata_file', "Gap_filled_meteorological_data_and_surface_energy.xml", ...
            'station_suffixes', "_surface.nc");
      case "firn_temperature"
         spec = struct( ...
            'product', product, ...
            'doi', "10.18739/A2833N00P", ...
            'metadata_file', "Firn_temperatures_and_measurement_depths_at_nine.xml", ...
            'station_suffixes', "_T_firn_obs.nc");
      case "simulated_firn"
         spec = struct( ...
            'product', product, ...
            'doi', "10.18739/A2CV4BS43", ...
            'metadata_file', "Simulated_firn_density_temperature_liquid_water.xml", ...
            'station_suffixes', ["_T_ice_bin_", "_rho_bin_", ...
            "_slwc_bin_", "_compaction_bin_"]);
   end

   % Fetch validation collects only NetCDF candidates through the shared fetch
   % helper. Surface/observed-firn names are complete filenames; simulated-firn
   % suffixes are prefixes followed by one numeric bin id and the .nc extension.
   spec.data_patterns = ["*.nc", fullfile("**", "*.nc")];
   if product == "simulated_firn"
      spec.station_file_mode = "numbered";
   else
      spec.station_file_mode = "exact";
   end
end
