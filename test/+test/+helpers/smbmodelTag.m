function tag = smbmodelTag(smbmodel)
%SMBMODELTAG Return canonical smbmodel tag for filenames and identifiers.
   arguments
      smbmodel = "all"
   end

   smbmodel = string(smbmodel);
   if any(strcmpi(smbmodel, "all"))
      tag = "all";
   else
      tag = test.helpers.sanitizeTag(join(smbmodel, "_"));
   end
end
