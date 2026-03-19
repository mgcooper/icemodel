function tf = isblankpart(value)
   %ISBLANKPART Return true for blank fullfile-style path parts.
   %
   %  TF = icemodel.isblankpart(VALUE) returns true when VALUE is blank text or
   %  an empty cell array. This is useful when composing optional path parts
   %  for ICEMODEL output and test-helper paths.
   %
   % See also icemodel.isblankinput isblanktext

   narginchk(1, 1)

   if iscell(value)
      tf = isempty(value);
   else
      tf = isblanktext(value);
   end
end
