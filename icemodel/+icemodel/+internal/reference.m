function str = reference(new)
   %REFERENCE set or get the icemodel reference
   %
   %  REF = ICEMODEL.INTERNAL.REFERENCE() Returns a char for the primary
   %  reference for icemodel.
   %
   %  REF = ICEMODEL.INTERNAL.REFERENCE(NEW_REFERENCE) Updates the primary
   %  reference for icemodel.
   %
   % See also: ICEMODEL.INTERNAL.VERSION

   % Return the icemodel reference or set a new one if supplied
   persistent current
   if isempty(current) || nargin == 1 && strcmp('reset', new)
      current = 'Cooper (2020), https://escholarship.org/uc/item/2fc6r410';
   elseif nargin == 1 && isrow(new) && ischar(new)
      current = new;
   end
   str = current;
end
