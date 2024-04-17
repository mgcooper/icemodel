function mkfolders(opts)

   simyears = opts.simyears;

   for MM = 1:numel(simyears)

      % set the output folder name and create it if it doesn't exist
      pathout = fullfile(opts.pathoutput, int2str(simyears(MM)));

      if ~isfolder(pathout)
         warning( ...
            'Output folders do not exist, creating them in %s', pathout)
         mkdir(pathout);
      end
   end
   
   % make a folder to save the model options
   if ~isfolder(fullfile(opts.pathoutput, 'opts'))
      mkdir(fullfile(opts.pathoutput, 'opts'));
   end
end
