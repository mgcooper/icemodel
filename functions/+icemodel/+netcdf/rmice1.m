function rmice1(pathdata, rmflag)

   arguments
      pathdata (1, :) char {mustBeFolder}
      rmflag (1, 1) logical {mustBeNumericOrLogical} = false
   end

   if rmflag == false
      return
   end

   allyears = 2009:2018;

   for n = 1:numel(allyears)

      thisyear = num2str(allyears(n));

      filelist = listfiles(fullfile(pathdata, thisyear), ...
         pattern="ice1_", aslist=true, fullpath=true);


   end
end
