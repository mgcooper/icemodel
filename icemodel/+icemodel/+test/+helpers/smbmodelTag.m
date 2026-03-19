function tag = smbmodelTag(smbmodel)
   %SMBMODELTAG Return canonical smbmodel tag for filenames and identifiers.
   arguments
      smbmodel = "all"
   end

   % Normalize mixed caller inputs before building the output tag.
   smbmodel = string(smbmodel);
   % Treat the aggregate selector specially; all other cases compose a tag.
   if any(strcmpi(smbmodel, "all"))
      tag = "all";
   else
      tag = icemodel.test.helpers.sanitizeTag(join(smbmodel, "_"));
   end
end
