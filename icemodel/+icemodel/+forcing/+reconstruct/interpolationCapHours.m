function caps = interpolationCapHours()
   %INTERPOLATIONCAPHOURS Approved per-channel interpolation ceilings.
   %
   %  caps = icemodel.forcing.reconstruct.interpolationCapHours()
   %
   % The default six-hour ceiling applies except where observed-only
   % holdouts support a channel rule: SWD and RH use nine hours, and
   % albedo uses 30 hours (D-39/D-42/D-50). D-49 applies SWD's same
   % nine-hour ceiling at a calendar-season boundary; the alias remains in
   % the returned contract so callers cannot drift to an independent value.

   caps = struct('default', 6, 'swd', 9, 'rh', 9, 'albedo', 30);
   caps.swd_season_boundary = caps.swd;
end
