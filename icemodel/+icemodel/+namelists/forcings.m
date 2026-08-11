function list = forcings()
   %FORCINGS Return the supported forcing-source names.
   %
   %  list = icemodel.namelists.forcings()
   %
   %  "promice" is the PROMICE AWS station-met source: a met_<site>_promice file
   %  that the verification staging tools build, where the SITENAME argument
   %  names the site. The per-station convention (forcings == sitename,
   %  e.g. "kanm" -> met_kanm_kanm) also works.
   %
   %  "promice_filled" is the CANONICAL runnable PROMICE forcing
   %  (met_<site>_promice_filled). The icemodel.forcing.reconstruct engine
   %  produces it with per-sample provenance. Native "promice" stays unmodified
   %  for provenance and QC. Its record is incomplete for most station-years, so
   %  it cannot force the model there.
   %
   %  "gcnet" is the separate Vandecrux gap-filled GC-Net surface/SEB source
   %  used for RetMIP Dye-2-long and Summit native forcing. "imau", "retmip",
   %  and "esm_snowmip" are native verification-staged runtime sources that the
   %  dataset importers write.

   list = ["mar"; "mar3.11"; "racmo"; "racmo2.3p3"; "merra"; "merra2"; ...
      "kanm"; "kanl"; "promice"; "promice_filled"; "gcnet"; "imau"; ...
      "retmip"; "esm_snowmip"];
end
