function opts = prepareRunOutput(opts)
   %PREPARERUNOUTPUT Prepare output folders and optional opts logging.
   %
   %  opts = icemodel.prepareRunOutput(opts)

   if ~opts.saveflag && ~(isfield(opts, 'saveopts') && opts.saveopts)
      return
   end

   if opts.saveflag
      % Create folders for each simulation year in output/ if they don't exist.
      icemodel.mkfolders(opts);
   elseif ~(exist(fullfile(opts.pathoutput, 'opts'), 'dir') == 7)
      mkdir(fullfile(opts.pathoutput, 'opts'));
   end

   if isfield(opts, 'saveopts') && opts.saveopts
      icemodel.saveRunOpts(opts);
   end
end
