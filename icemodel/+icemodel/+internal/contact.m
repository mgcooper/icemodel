function str = contact(new)
   %CONTACT set or get the icemodel contact
   %
   %  STR = ICEMODEL.INTERNAL.CONTACT() Returns a char for the primary
   %  contact for icemodel.
   %
   %  STR = ICEMODEL.INTERNAL.CONTACT(NEW_CONTACT) Updates the primary
   %  contact for icemodel.
   %
   % See also: ICEMODEL.INTERNAL.VERSION

   persistent current
   if isempty(current) || nargin == 1 && strcmp('reset', new)
      current = 'matt dot cooper at pnnl dot gov';
   elseif nargin == 1 && isrow(new) && ischar(new)
      current = new;
   end
   str = current;
end
