function list = userdata()
   %USERDATA Return the supported core userdata source names.
   %
   %  list = icemodel.namelists.userdata()

   % The list includes promice_filled because setopts defaults
   % userdata = forcings. The same-name guard in loadmet then disables
   % swapping, which is the intended no-op for the gap-filled product.
   list = ["mar"; "mar3.11"; "modis"; "merra"; "merra2"; ...
      "racmo"; "racmo2.3p3"; "kanm"; "kanl"; "promice"; ...
      "promice_filled"; "gcnet"; "imau"; "retmip"; "esm_snowmip"];
end
