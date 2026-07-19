function albedo = normalizeGeusModisAlbedo(albedo)
   %NORMALIZEGEUSMODISALBEDO Mask undocumented GEUS albedo sentinels.
   %
   %  albedo = icemodel.forcing.helpers.normalizeGeusModisAlbedo(albedo)
   %
   % GEUS Greenland Reflectivity C6 files use finite 999 values for missing
   % albedo without declaring a NetCDF fill or valid-range attribute. Convert
   % every nonfinite or out-of-domain sample to NaN at the shared source boundary
   % so point, polygon, builder, and bounded-repair paths agree on coverage.

   albedo = double(albedo);
   invalid = ~isfinite(albedo) | imag(albedo) ~= 0 ...
      | real(albedo) < 0 | real(albedo) > 1;
   albedo(invalid) = NaN;
end
