function tf = istextlike(x, varargin)
   %ISTEXTLIKE Return true when input contains only text values.
   %
   %  TF = ISTEXTLIKE(X) returns true when X is text-like:
   %   - a char row vector, including ''
   %   - a string array
   %   - a cell array whose elements are themselves text-like
   %
   %  TF = ISTEXTLIKE(X, MODE) controls how cell contents are reduced. MODE may
   %  be 'all', 'any', or 'each'. The default is 'all'.
   %
   %  TF = ISTEXTLIKE(..., 'nontrivial') requires nonzero-length text values.
   %
   %  This is the container-level text predicate used by icemodel's
   %  dependency-layer validators. Use ISSCALARTEXT for legacy scalar parsing
   %  semantics and ISBLANKTEXT for blank-text detection.
   %
   % See also ischarlike isscalartext isblanktext
   %
   %#codegen

   [opt, args, nargs] = parseoptarg(varargin, {'all', 'any', 'each'}, 'all');

   if nargs > 0 && strcmpi(args{1}, 'nontrivial')
      hastext = @(value) all(strlength(string(value)) > 0, 'all');
   else
      hastext = @(value) true;
   end

   if iscell(x)
      switch opt
         case 'all'
            tf = all(cellfun(@(arg) istextlike(arg, 'all', args{:}), x));
         case 'any'
            tf = any(cellfun(@(arg) istextlike(arg, 'all', args{:}), x));
         case 'each'
            tf = cellfun(@(arg) istextlike(arg, 'all', args{:}), x);
      end
      return
   end

   if isstring(x)
      tf = hastext(x);
   else
      tf = ischar(x) && (isrow(x) || isequal(x, '')) && hastext(x);
   end
end
