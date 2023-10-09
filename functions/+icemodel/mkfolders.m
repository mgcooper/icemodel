function mkfolders(opts)

simyears = opts.simyears;

for MM = 1:numel(simyears)
   
   % set the output folder name and create it if it doesn't exist
   pathout = fullfile(opts.pathoutput, int2str(simyears(MM)));
   
   if ~isfolder(pathout)
      mkdir(pathout); 
   end
   
end