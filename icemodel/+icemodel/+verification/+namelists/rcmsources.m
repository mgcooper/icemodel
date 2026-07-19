function list = rcmsources()
   %RCMSOURCES Verification RCM source labels in canonical staging order.
   %
   %  list = icemodel.verification.namelists.rcmsources()
   %
   % This list is shared by importers, RCM staging helpers, and manifest
   % source-list derivation so the model set and order cannot drift.
   list = ["mar", "merra", "racmo"];
end
