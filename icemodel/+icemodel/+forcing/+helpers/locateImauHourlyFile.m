function filename = locateImauHourlyFile(source_dir, station)
   %LOCATEIMAUHOURLYFILE Find a flat or hourly-subfolder PANGAEA tab file.
   patterns = [ ...
      fullfile(source_dir, 'hourly', "*" + station + "*.tab")
      fullfile(source_dir, "*" + station + "*.tab")
      fullfile(source_dir, '**', "*" + station + "*.tab")];
   for pattern = reshape(patterns, 1, [])
      hits = dir(pattern);
      hits = hits(~[hits.isdir]);
      if isempty(hits)
         continue
      end
      folders = string({hits.folder});
      names = string({hits.name});
      keep = ~contains(folders, filesep + "daily") ...
         & ~startsWith(names, "GRL_" + station + "_AWS");
      hits = hits(keep);
      if ~isempty(hits)
         filename = string(fullfile(hits(1).folder, hits(1).name));
         return
      end
   end
   error('icemodel:forcing:locateImauHourlyFile:fileNotFound', ...
      'IMAU hourly source file for %s not found under %s', station, source_dir)
end
