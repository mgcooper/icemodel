function pathlist = setpath(pathtype, sitename, simmodel, userdata, simyears, varargin)
   %SETPATH sets icemodel input and output paths
   %
   %

   narginchk(2, Inf)

   if nargin < 5 || isempty(simyears)
      simyears = '';
   else
      simyears = arrayfun(@num2str, simyears(:), 'un', false);
   end

   if nargin < 4, userdata = ''; end
   if nargin < 3, simmodel = ''; end

   switch pathtype
      case 'output'
         pathlist = fullfile(getenv('ICEMODELOUTPUTPATH'), sitename, simmodel, ...
            userdata, simyears, varargin{:});

      case 'input'
   end
end