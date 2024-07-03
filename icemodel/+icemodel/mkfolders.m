function mkfolders(opts)

   simyears = opts.simyears;

   for MM = 1:numel(simyears)

      % set the output folder name and create it if it doesn't exist
      pathout = fullfile(opts.pathoutput, int2str(simyears(MM)));

      if ~(exist(pathout, 'dir') == 7)
         warning( ...
            'Output folders do not exist, creating them in %s', pathout)
         mkdir(pathout);
      end
   end

   % make a folder to save the model options
   if ~(exist(fullfile(opts.pathoutput, 'opts'), 'dir') == 7)
      mkdir(fullfile(opts.pathoutput, 'opts'));
   end
end
