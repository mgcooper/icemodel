function tf = isblanktext(x)
   %ISBLANKTEXT Return true for omitted or zero-length text inputs
   %
   %  TF = ISBLANKTEXT(X) returns TF = true when X is:
   %   1. an empty array such as string.empty()
   %   2. a scalar text value with zero characters, such as '' or ""
   %
   % This complements ISSCALARTEXT by distinguishing blank/omitted text from
   % nonblank scalar text.
   %
   % Examples
   % tf = isblanktext("")
   % tf =
   % logical
   %  1
   %
   % tf = isblanktext(string.empty())
   % tf =
   % logical
   %  1
   %
   % tf = isblanktext("abc")
   % tf =
   % logical
   %  0

   narginchk(1, 1)

   if isempty(x)
      tf = true;
      return
   end

   tf = isscalartext(x) && strlength(string(x)) == 0;
end
