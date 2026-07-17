function ids = rcmProductIds(labels)
   %RCMPRODUCTIDS Map RCM runtime/storage labels to explicit product ids.
   %
   %  ids = icemodel.verification.namelists.rcmProductIds(["mar","merra"])
   %
   % The short labels remain user-facing model-family selectors. Staged
   % met/userdata artifacts and manifest source lists use these explicit product
   % ids so future versions from the same RCM family can coexist.

   arguments
      labels (1, :) string = icemodel.verification.namelists.rcmsources()
   end

   ids = strings(size(labels));
   for k = 1:numel(labels)
      switch labels(k)
         case "mar"
            ids(k) = "mar3.11";
         case "merra"
            ids(k) = "merra2";
         case "racmo"
            ids(k) = "racmo2.3p3";
         otherwise
            ids(k) = labels(k);
      end
   end
end
