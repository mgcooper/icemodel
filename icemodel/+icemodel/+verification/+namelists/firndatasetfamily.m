function families = firndatasetfamily()
   %FIRNDATASETFAMILY Return verification families used by firn staging previews.
   %
   %  families = icemodel.verification.namelists.firndatasetfamily()
   %
   % These families participate in the firn model-development staging workflow.
   % Snow-only ESM-SnowMIP and Laugh-test families remain in datasetfamily().

   families = ["promice", "retmip", "imau", "research_site", "sumup"];
end
