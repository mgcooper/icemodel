function list = userdata()
%USERDATA Return the supported core userdata source names.
%
%  list = icemodel.namelists.userdata()

   list = ["mar"; "mar3.11"; "modis"; "merra"; "merra2"; ...
      "racmo"; "racmo2.3p3"; "kanm"; "kanl"; "promice"; ...
      "gcnet"; "imau"; "retmip"; "esm_snowmip"];
end
