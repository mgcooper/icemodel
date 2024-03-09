function varargout = getdefaults(whichfile, whichprops, output_type)

   arguments
      whichfile (1, :) char {mustBeMember(whichfile, ...
         {'ice1', 'ice2', 'met', 'grid'})} ...
         = 'ice1'
      whichprops (1, :) string {mustBeMember(whichprops, ...
         ["varnames", "longnames", "units"])} ...
         = ["varnames", "longnames", "units"]
      output_type (1, :) char {mustBeMember(output_type, ...
         {'asstruct', 'aslist'})} ...
         = 'asstruct'
   end
   
   % Overrule output_type if nargout matches the number of requested props
   if nargout == numel(whichprops)
      output_type = 'aslist';
   end

   for thisprop = whichprops(:)'
      propvals.(thisprop) = icemodel.netcdf.defaults.(thisprop)(whichfile);
   end
   
   if output_type == "aslist"
      nargoutchk(numel(whichprops), numel(whichprops));
      varargout = struct2cell(propvals);
      
   else
      nargoutchk(1, 1);
      varargout{1} = propvals;
   end
