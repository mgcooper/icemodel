function removeBaselineProfileStage(pathname)
   %REMOVEBASELINEPROFILESTAGE Remove one owned profiler staging directory.
   %
   %  icemodel.test.helpers.removeBaselineProfileStage(pathname)

   arguments
      pathname string {mustBeTextScalarOrEmpty}
   end

   if ~isblanktext(pathname) && isfolder(pathname)
      rmdir(pathname, 's');
   end
end
