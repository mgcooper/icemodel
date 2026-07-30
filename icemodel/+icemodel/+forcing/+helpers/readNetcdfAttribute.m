function value = readNetcdfAttribute(filename, variable, attribute)
   %READNETCDFATTRIBUTE Return one NetCDF attribute or empty string when absent.
   try
      value = string(ncreadatt(filename, variable, attribute));
   catch
      value = "";
   end
end
