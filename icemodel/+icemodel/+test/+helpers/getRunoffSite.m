function runoff_site = getRunoffSite(sitename)
   %GETRUNOFFSITE Map formal station cases to runoff-validation catchments.

   % Normalize the incoming selector so the mapping works for any text type.
   sitename = string(sitename);
   runoff_site = sitename;

   % Preserve direct names unless the formal station maps to a catchment.
   runoff_site(lower(sitename) == "kanm") = "behar";
   runoff_site(lower(sitename) == "kanl") = "upperbasin";
end
