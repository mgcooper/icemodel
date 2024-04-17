function varargout = getdefaults(whichfiles, whichprops, output_type)

   arguments
      whichfiles (1, :) string {mustBeMember(whichfiles, ...
         {'ice1', 'ice2', 'met', 'dims', 'dimensions'})} ...
         = 'ice1'
      whichprops (1, :) string {mustBeMember(whichprops, ...
         ["varnames", "standardnames", "longnames", "units", "axes"])} ...
         = ["varnames", "standardnames", "longnames", "units", "axes"]
      output_type (1, :) char {mustBeMember(output_type, ...
         {'asstruct', 'aslist', 'astable'})} ...
         = 'asstruct'
   end

   % Overrule output_type if nargout matches the number of requested props
   if nargout == numel(whichprops)
      output_type = 'aslist';
   end

   for thisfile = whichfiles(:)'

      for thisprop = whichprops(:)'
         props.(thisprop).(thisfile) = ...
            icemodel.netcdf.defaults.(thisprop)(thisfile);
      end

      % Remove 'depth' dim from ice1 - this is done for defdimvars
      props = dropdims(props, whichprops, thisfile);
   end

   varargout = struct2cell(props);

   % if output_type == "aslist"
   %    nargoutchk(numel(whichprops), numel(whichprops));
   %    varargout = struct2cell(props);
   % end
end

%% drop dims
function props = dropdims(props, propnames, datafile)

   if strcmp(datafile, 'ice1')
      drop = ismember(props.varnames.dims, 'depth');
      for thisprop = propnames(:)'
         props.(thisprop).dims(drop) = [];
      end
   end
end
% function [varnames, axes, units, longnames, standardnames] = dropdims( ...
%       datafile, varnames, axes, units, longnames, standardnames)
%
%    if strcmp(datafile, 'ice1')
%       drop = ismember(varnames.dims, 'depth');
%       varnames.dims(drop) = [];
%       axes.dims(drop) = [];
%       units.dims(drop) = [];
%       longnames.dims(drop) = [];
%       standardnames.dims(drop) = [];
%    end
% end
