function casename = setcase(forcings, userdata, uservars, ...
      smbmodel, sitename, simyear)

   % Input parsing
   if nargin < 4 || isempty(smbmodel); smbmodel = ''; end
   if nargin < 5 || isempty(sitename); sitename = ''; end
   if nargin < 6 || isempty(simyear); simyear = ''; end

   if isnumeric(simyear)
      simyear = num2str(simyear);
   end

   % convertStringsToChars in a pre-R2017b compatible way:
   args = {forcings, userdata, uservars, smbmodel, sitename, simyear};
   for n = 1:numel(args)
      if isstring(args{n})
         args{n} = char(args{n});
      end
   end
   [smbmodel, sitename, forcings, userdata, uservars, simyear] = deal(args{:});

   % Set the casename
   if nargin == 3
      casename = [forcings '_forcings_' upper(userdata) '_' uservars];

   elseif nargin == 6

      % The legacy filename format:
      casename = [smbmodel '_' sitename '_' simyear '_' upper(forcings) ...
         '_forcings_' upper(userdata) '_' uservars '.mat'];
   end
end
