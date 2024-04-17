function casename = setcase(forcings, userdata, uservars, ...
      simmodel, sitename, simyear)
   arguments
      forcings char
      userdata char
      uservars char
      simmodel char = ''
      sitename char = ''
      simyear char = ''
   end
   if nargin == 3
      casename = [forcings '_swap_' upper(userdata) '_' uservars];

   elseif nargin == 6

      % The legacy filename format:
      if isnumeric(simyear)
         simyear = num2str(simyear);
      end

      casename = [simmodel '_' sitename '_' simyear '_' upper(forcings) ...
         '_swap_' upper(userdata) '_' uservars '.mat'];
   end
end
