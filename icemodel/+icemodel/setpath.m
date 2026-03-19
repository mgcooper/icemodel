function pathlist = setpath(varargin)
   %SETPATH Compatibility wrapper for icemodel.getpath.
   %
   %  pathlist = icemodel.setpath(...)
   %
   % Use icemodel.getpath(...) for new code. This wrapper exists so older
   % callers continue to resolve the same canonical paths while the repo
   % transitions to the clearer getter name.

   pathlist = icemodel.getpath(varargin{:});
end
