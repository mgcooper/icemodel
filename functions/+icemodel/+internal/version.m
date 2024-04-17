function str = version(new)
   %VERSION set or get the icemodel version number
   %
   %  VERS = ICEMODEL.INTERNAL.VERSION()
   %  VERS = ICEMODEL.INTERNAL.VERSION(NEW_VERSION)
   %
   % See also: ICEMODEL.INTERNAL.REFERENCE

   % Return the icemodel version or set a new one if supplied
   persistent current
   if isempty(current) || nargin == 1 && strcmp('reset', new)
      current = '1.0.0';
   elseif nargin == 1 && isrow(new) && ischar(new)
      current = new;
   end
   str = current;
end
