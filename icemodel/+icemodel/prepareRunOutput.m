function opts = prepareRunOutput(opts)
   %PREPARERUNOUTPUT Prepare output folders and optional opts logging.
   %
   %  opts = icemodel.prepareRunOutput(opts)

   if ~opts.saveflag ...
         && ~(isfield(opts, 'saveopts') && opts.saveopts) ...
         && ~(isfield(opts, 'saverestart') && opts.saverestart)
      return
   end

   if opts.saveflag
      % Create folders for each simulation year in output/ if they don't exist.
      icemodel.mkfolders(opts);
   elseif isfield(opts, 'saveopts') && opts.saveopts ...
         && ~isfolder(fullfile(opts.pathoutput, 'opts'))
      mkdir(fullfile(opts.pathoutput, 'opts'));
   end

   if isfield(opts, 'saveopts') && opts.saveopts
      icemodel.saveRunOpts(opts);
   end

   if isfield(opts, 'saverestart') && opts.saverestart ...
         && ~isfolder(opts.pathrestart)
      mkdir(opts.pathrestart);
   end
end
