function list = tier()
   %TIER Return the supported snow-verification suite tiers.
   %
   %  list = icemodel.verification.namelists.tier()
   %
   % Outputs
   %  list   String column of supported suite tier selectors.
   %
   % Role
   %  Canonical selector list shared by setup importers, validators, and
   %  normal verification workflow filters.

   % "full" is retained as a selector even before full-size fixtures are staged.
   list = ["smoke"; "full"];
end
