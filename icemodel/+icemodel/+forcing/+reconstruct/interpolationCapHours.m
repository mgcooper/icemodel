function caps = interpolationCapHours()
   %INTERPOLATIONCAPHOURS Approved per-channel interpolation ceilings.
   %
   %  caps = icemodel.forcing.reconstruct.interpolationCapHours()
   %
   % The default six-hour ceiling applies except where observed-only
   % holdouts support a channel rule: SWD and RH use nine hours, and
   % albedo uses 30 hours (D-39/D-42/D-50). D-49 applies the same nine-hour
   % SWD ceiling at a calendar-season boundary. The returned struct carries
   % that alias, so every caller reads one value.

   caps = struct('default', 6, 'swd', 9, 'rh', 9, 'albedo', 30);
   caps.swd_season_boundary = caps.swd;
end
