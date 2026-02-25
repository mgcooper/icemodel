function varargout = getopts(opts, varargin)
   %GETOPTS Return model options by name.
   %
   % [a, b] = icemodel.getopts(opts, 'a', 'b')
   % Input option names define return order.
   %
   %#codegen

   narginchk(1, inf);

   if (nargin == 1 && nargout == 1) ...
         || (nargin > 1 && strcmp('all', varargin{1}))
      varargout{1} = opts;
      return
   end

   for n = 1:nargin-1
      arg = varargin{n};
      varargout{n} = opts.(arg);
   end
end
