function pathlist = setpath(varargin)
   %SETPATH Compatibility wrapper for icemodel.getpath.
   %
   %  pathlist = icemodel.setpath(...)
   %
   % Use icemodel.getpath(...) for new code. This wrapper resolves the same
   % canonical paths for callers that still use the setpath name.

   pathlist = icemodel.getpath(varargin{:});
end
