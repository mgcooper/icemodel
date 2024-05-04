function casename = setcase(forcings, userdata, uservars, ...
      smbmodel, sitename, simyear)
   arguments
      forcings char
      userdata char
      uservars char
      smbmodel char = ''
      sitename char = ''
      simyear char = ''
   end
   if nargin == 3
      casename = [forcings '_forcings_' upper(userdata) '_' uservars];

   elseif nargin == 6

      % The legacy filename format:
      if isnumeric(simyear)
         simyear = num2str(simyear);
      end

      casename = [smbmodel '_' sitename '_' simyear '_' upper(forcings) ...
         '_forcings_' upper(userdata) '_' uservars '.mat'];
   end
end
