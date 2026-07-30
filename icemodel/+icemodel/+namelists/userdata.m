function list = userdata()
%USERDATA Return the supported core userdata source names.
%
%  list = icemodel.namelists.userdata()

   % promice_filled is listed because setopts defaults userdata = forcings;
   % loadmet's same-name guard then disables swapping, which is the intended
   % no-op for the gap-filled product.
   list = ["mar"; "mar3.11"; "modis"; "merra"; "merra2"; ...
      "racmo"; "racmo2.3p3"; "kanm"; "kanl"; "promice"; ...
      "promice_filled"; "gcnet"; "imau"; "retmip"; "esm_snowmip"];
end
