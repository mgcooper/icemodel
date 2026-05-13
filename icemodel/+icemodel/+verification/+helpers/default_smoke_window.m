function [window_start, window_end] = default_smoke_window(sitename)
   %DEFAULT_SMOKE_WINDOW Per-site default verification window (one snow water year).
   %
   %  [start, end] = icemodel.verification.helpers.default_smoke_window()
   %  [start, end] = icemodel.verification.helpers.default_smoke_window("cdp")
   %
   %  Returns the canonical "one snow water year" window for the given
   %  ESM-SnowMIP site, used as the default staging window for
   %  importEsmSnowmip and the default runtime window for
   %  run_snow_verification_suite. The window is the second insitu year
   %  (start_year + 1) covering [Oct 1, 0:00 UTC] to [Sep 30, 23:00 UTC].
   %  Picking the second year lets staging skip any leading spin-up gaps
   %  in the upstream PANGAEA files.
   %
   %  With no argument, defaults to "cdp" (Col de Porte, Menard 2019
   %  ESSD) — the most canonical / widely-cited ESM-SnowMIP snow
   %  verification site. Single-site single-year default keeps
   %  interactive verification runs fast.
   %
   %  Inputs
   %    sitename : string  ESM-SnowMIP site code. Default "cdp".
   %
   %  Returns
   %    window_start : datetime (TimeZone='UTC')  Oct 1 00:00 of the
   %                                              second insitu year.
   %    window_end   : datetime (TimeZone='UTC')  Sep 30 23:00 of the
   %                                              following year.
   %
   % See also: icemodel.verification.helpers.snowmipinfo,
   %  icemodel.verification.setup.importEsmSnowmip,
   %  run_snow_verification_suite

   arguments
      sitename (1, 1) string ...
         {icemodel.verification.validators.mustBeSnowmipSite} = "cdp"
   end

   info = icemodel.verification.helpers.snowmipinfo(sitename);
   smoke_year = info.insitu_window(1) + 1;
   window_start = datetime(smoke_year, ...
      10, 1, 0, 0, 0, 'TimeZone', 'UTC');
   window_end = datetime(smoke_year + 1, ...
      9, 30, 23, 0, 0, 'TimeZone', 'UTC');
end
