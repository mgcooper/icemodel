function list = forcings()
%FORCINGS Return the supported forcing-source names.
%
%  list = icemodel.namelists.forcings()
%
%  "promice" is the generic PROMICE/GC-Net AWS station-met source: a
%  met_<site>_promice file built by the verification staging tools, where the
%  site is named by the SITENAME argument. The legacy per-station convention
%  (forcings == sitename, e.g. "kanm" -> met_kanm_kanm) is still supported.

   list = ["mar"; "racmo"; "merra"; "kanm"; "kanl"; "promice"];
end
