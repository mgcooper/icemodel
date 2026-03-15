function runoff_site = getRunoffSite(sitename)
%GETRUNOFFSITE Map formal station cases to runoff-validation catchments.

   sitename = string(sitename);
   runoff_site = sitename;

   runoff_site(lower(sitename) == "kanm") = "behar";
   runoff_site(lower(sitename) == "kanl") = "upperbasin";
end
