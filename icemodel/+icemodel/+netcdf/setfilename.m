function filename = setfilename(datafile, smbmodel, forcings, userdata, ...
      sitename, simyears, filepath, opts)

   arguments
      datafile (1, :) char {mustBeMember(datafile, ...
         {'ice1', 'ice2', 'met'})}
      smbmodel (1, :) char {mustBeTextScalar}
      forcings (1, :) char {mustBeTextScalar}
      userdata (1, :) char {mustBeTextScalar}
      sitename (1, :) char {mustBeTextScalar}
      simyears (1, :) string % casts double to string
      filepath (1, :) char {mustBeFolder} = getenv('ICEMODELOUTPUTPATH')
      opts.extension (1, :) char {mustBeMember(opts.extension, {'nc', 'nc4'})} = 'nc4'
   end

   filename = cell(numel(simyears), 1);

   for n = 1:numel(simyears)

      filename{n} = fullfile(filepath, strjoin( ...
         {smbmodel, datafile, forcings, userdata, sitename, char(simyears(n)), ...
         opts.extension}, '.'));
   end

   if numel(simyears) == 1
      filename = filename{1};
   end

   % temporary hack - replace 'sector' with 'sw'
   filename = strrep(filename, '.sector.', '.sw.');
end
