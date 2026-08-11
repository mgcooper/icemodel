function list = rcmsources()
   %RCMSOURCES Verification RCM source labels in canonical staging order.
   %
   %  list = icemodel.verification.namelists.rcmsources()
   %
   % Importers, RCM staging helpers, and manifest source-list derivation all
   % read this list, so they use the same model set in the same order.
   list = ["mar", "merra", "racmo"];
end
