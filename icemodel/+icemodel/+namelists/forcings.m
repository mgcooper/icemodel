function list = forcings()
%FORCINGS Return the supported forcing-source names.
%
%  list = icemodel.namelists.forcings()
%
%  "promice" is the PROMICE AWS station-met source: a met_<site>_promice file
%  built by the verification staging tools, where the site is named by the
%  SITENAME argument. The legacy per-station convention (forcings == sitename,
%  e.g. "kanm" -> met_kanm_kanm) is still supported. "promice_filled" is the
%  CANONICAL runnable PROMICE forcing (met_<site>_promice_filled), produced
%  by the icemodel.forcing.reconstruct engine with per-sample provenance;
%  native "promice" is retained unmodified for provenance and QC — its
%  record is incomplete for most station-years and cannot force the model
%  there. "gcnet" is the distinct Vandecrux gap-filled GC-Net surface/SEB source
%  used for RetMIP Dye-2-long and Summit native forcing. "imau", "retmip",
%  and "esm_snowmip" are native verification-staged runtime sources written
%  by the dataset importers.

   list = ["mar"; "mar3.11"; "racmo"; "racmo2.3p3"; "merra"; "merra2"; ...
      "kanm"; "kanl"; "promice"; "promice_filled"; "gcnet"; "imau"; ...
      "retmip"; "esm_snowmip"];
end
