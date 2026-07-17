function list = forcings()
%FORCINGS Return the supported forcing-source names.
%
%  list = icemodel.namelists.forcings()
%
%  "promice" is the PROMICE AWS station-met source: a met_<site>_promice file
%  built by the verification staging tools, where the site is named by the
%  SITENAME argument. The legacy per-station convention (forcings == sitename,
%  e.g. "kanm" -> met_kanm_kanm) is still supported. "gcnet" is the distinct
%  Vandecrux gap-filled GC-Net surface/SEB source used for RetMIP Dye-2-long and
%  Summit native forcing. "imau", "retmip", and "esm_snowmip" are native
%  verification-staged runtime sources written by the dataset importers.

   list = ["mar"; "mar3.11"; "racmo"; "racmo2.3p3"; "merra"; "merra2"; ...
      "kanm"; "kanl"; "promice"; "gcnet"; "imau"; "retmip"; ...
      "esm_snowmip"];
end
