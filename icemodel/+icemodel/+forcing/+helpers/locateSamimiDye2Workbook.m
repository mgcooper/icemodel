function filename = locateSamimiDye2Workbook(source_dir)
   %LOCATESAMIMIDYE2WORKBOOK Find the Dye-2 summer 2016 AWS workbook.
   preferred = fullfile(source_dir, "Dye2_AWS_Summer2016.xlsx");
   if isfile(preferred)
      filename = string(preferred);
      return
   end
   files = icemodel.forcing.helpers.listSourceFiles(source_dir);
   [~, names, extensions] = fileparts(files);
   basenames = names + extensions;
   keep = contains(basenames, "Dye", 'IgnoreCase', true) ...
      & contains(basenames, "AWS", 'IgnoreCase', true) ...
      & contains(basenames, "2016", 'IgnoreCase', true) ...
      & endsWith(basenames, ".xlsx", 'IgnoreCase', true);
   hits = files(keep);
   if isempty(hits)
      error('icemodel:forcing:locateSamimiDye2Workbook:fileNotFound', ...
         'Samimi Dye-2 AWS workbook not found under %s', source_dir)
   end
   filename = string(hits(1));
end
