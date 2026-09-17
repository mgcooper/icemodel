function tf = islogicalscalar(x)
   %ISLOGICALSCALAR Determine if input is a scalar logical value.
   %
   %  TF = ISLOGICALSCALAR(X) returns true if X is a logical value with
   %  exactly one element, and false otherwise. A logical array with more
   %  than one element returns false, as do an empty logical, a numeric 0
   %  or 1, and a char.
   %
   % Example
   %  islogicalscalar(true)          % 1
   %  islogicalscalar([true false])  % 0
   %  islogicalscalar(1)             % 0
   %
   % See also: islogical, isscalar, isscalartext
   %
   %#codegen

   % parseoptarg calls this to tell a logical DEFAULTOPT from a text one.
   % Both must hold: logical type and exactly one element.
   tf = islogical(x) && isscalar(x);
end
