function out = tocellstr(in, varargin)
   %TOCELLSTR Convert char, string arrays, or cell arrays to cellstr.
   %
   %  OUT = TOCELLSTR(IN)
   %  OUT = TOCELLSTR(IN, DATEFMT)
   %
   % Description
   %  OUT = TOCELLSTR(IN) tries to convert IN to a cellstr array. If it fails,
   %  it returns IN unchanged.
   %
   %  OUT = TOCELLSTR(IN, DATEFMT) additionally supplies DATEFMT to cellstr to
   %  convert datetime or duration IN to date formatted text. See `cellstr`
   %  documentation for details.
   %
   % See also: cellstr tocolumn torow rowvec
   %
   %#codegen

   try
      out = cellstr(in, varargin{:});
   catch e
      wid = ['custom:' mfilename ':conversionToCellstrFailed'];
      msg = 'Conversion to cellstr failed. Returning input unchanged.';
      warning(wid, msg)
      out = in;
   end
end
