function tf = isblanktext(x)
   %ISBLANKTEXT Return true for blank text values.
   %
   %  TF = ISBLANKTEXT(X) returns TF = true when X is blank text:
   %   1. a zero-length char array such as '' or char.empty()
   %   2. a zero-length string value such as "" or string.empty()
   %
   %  TF = ISBLANKTEXT(X) returns false for non-text inputs, missing strings,
   %  nonscalar string arrays, and nonblank text values.
   %
   %  This complements ISSCALARTEXT by detecting blank text directly rather
   %  than inheriting ISSCALARTEXT's legacy parser behavior.
   %
   % See also isscalartext

   if ischar(x)
      tf = isempty(x) || (isrow(x) && strlength(string(x)) == 0);
      return
   end

   if isstring(x)
      tf = isempty(x) || (isscalar(x) && ~ismissing(x) && strlength(x) == 0);
      return
   end

   tf = false;
end
