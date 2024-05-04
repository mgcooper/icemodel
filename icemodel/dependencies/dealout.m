function varargout = dealout(varargin)
   %DEALOUT deal out function outputs
   %
   % [argout1, argout2] = dealout(argin1, argin2, ..., arginN]
   % [cellArrayOutput{1:nargout}] = dealout(cellArrayInput{:});
   % [-] = dealout(structInput);
   % [-] = dealout(struct2cell(structInput));
   %
   % The number of output arguments does not have to match the number of inputs,
   % but they will be dealt out in the exact order they are dealt in i.e.
   % argout1 = argin1, argout2 = argin2, and so forth.
   %
   % Also the number of outputs must of course not exceed the number of inputs.
   %
   % Note: if the calling function tries to do something like this:
   % varargout = dealout(arg1) ... and no output is requested from the base
   % workspace, then the first element of arg1 will get sent back. For example:
   %
   % function celloutput = myfunction(varargin)
   %
   % ... function code
   %
   % varargout = dealout(cellarg1);
   %
   % % Then, in a script:
   % myfunction(...)
   %
   % % with no requested outputs, will return the first element of cellarg1. But
   % if this syntax is used, nothing will be returned:
   %
   % function celloutput = myfunction(varargin)
   %
   % ... function code
   %
   % [varargout{1:nargout}] = dealout(cellarg1);
   %
   % Matt Cooper, 22 Jun 2023
   %
   % See also

   args = varargin;

   % This first part is designed to eliminate the calling syntax:
   %
   % [val1, val2, ..., valN] = dealout(struct2cell(opts))
   %
   % Here, opts is a name-value struct with N pairs. If that calling syntax is
   % used, struct2cell produces a cell-array with one element per name-value
   % pair in opts, which is the desired result, but upon passing it to this
   % function, it gets nested inside varargin, and then ARGS is a scalar
   % cell-array, with the desired cell-array its only element.
   %
   % In contrast, if struct2cell is applied here, as it is below, it produces
   % the desired cell-array and is not nested inside the scalar varargin array.
   %
   % Note that, if args = varargin{:} is used instead of args = varargin, it
   % would produce the desired result, but it is unclear what the side effects
   % of this would be in other use cases.

   if numel(args) == 1 && isstruct(args{1})
      try
         args = struct2cell(args{:});
      catch
      end
   end

   if nargout > numel(args) && nargout > numel(args{:})
      error('One input required for each requested output')
   end

   try
      % Syntax is [out1, out2, ..., outN] = dealout(in1, in2, ..., inN]
      [varargout{1:nargout}] = deal(args{1:nargout});

      % Note: this syntax is useful for forcing one output
      % [varargout{1:max(1,nargout)}] = deal(varargin{1:nargout});
   catch
      % Syntax is [out1, out2, ..., outN] = dealout(CellArray)
      varargout = args{:};
   end
end
