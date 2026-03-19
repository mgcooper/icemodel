function restart = loadRestartState(opts)
%LOADRESTARTSTATE Load a saved year-boundary restart state.
%
%  restart = icemodel.loadRestartState(opts)
%
% Requires:
%  opts.use_restart = true
%  opts.restartfile = full path to a saved restart state

   if ~isfield(opts, 'restartfile') || isempty(opts.restartfile) ...
         || isblanktext(opts.restartfile)
      error('opts.restartfile must be set when opts.use_restart is true')
   end

   filepath = char(opts.restartfile);
   if exist(filepath, 'file') ~= 2
      error('restart file does not exist: %s', filepath)
   end

   S = load(filepath, 'restart');
   if ~isfield(S, 'restart')
      error('restart file does not contain struct "restart": %s', filepath)
   end

   restart = S.restart;
end
