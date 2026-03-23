function list = uservars()
   %USERVARS Return the supported core userdata variable names.
   %
   %  list = icemodel.namelists.uservars()
   %
   % TODO: add all forcings e.g. ppt

   list = ["albedo"; "tair"; "swd"; "lwd"; "rh"; "wsdp"; "psfc"];
end
