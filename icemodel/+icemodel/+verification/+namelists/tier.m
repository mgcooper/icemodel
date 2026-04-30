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

   % "smoke" stages a single representative snow year per site (the
   % default agent-feedback loop), "full" stages the full available
   % insitu range, and "custom" requires explicit startdate / enddate
   % when staging via importEsmSnowmip.
   list = ["smoke"; "full"; "custom"];
end
