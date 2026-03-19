function tf = isblankinput(value)
   %ISBLANKINPUT Return true for omitted optional inputs.
   %
   %  TF = icemodel.isblankinput(VALUE) returns true when VALUE is either:
   %   1. an empty array such as []
   %   2. blank text such as '', "", or string.empty()
   %
   %  This helper preserves the older parser behavior used by optional inputs
   %  that historically treated [] as an omitted placeholder. Use ISBLANKTEXT
   %  when the contract is specifically about text values.
   %
   % See also isblanktext icemodel.isblankpart

   narginchk(1, 1)
   tf = isempty(value) || isblanktext(value);
end
