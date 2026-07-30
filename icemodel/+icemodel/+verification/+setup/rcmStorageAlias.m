function alias = rcmStorageAlias(dataset_family, case_id)
   %RCMSTORAGEALIAS Return the collision-safe RCM artifact identity for a case.
   %
   % Manifest case ids remain user-facing scientific identities. Only the two
   % proved cache collisions need distinct on-disk names; every other case keeps
   % its existing id exactly.

   arguments
      dataset_family (1, 1) string
      case_id (1, 1) string
   end

   alias = case_id;
   if dataset_family == "retmip" && lower(case_id) == "kanu"
      alias = "retmip_kanu";
   elseif dataset_family == "sumup" ...
         && lower(case_id) == "firn aquifer (fa)"
      alias = "fa";
   end
end
