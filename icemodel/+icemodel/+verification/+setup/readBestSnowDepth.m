function [values, source] = readBestSnowDepth(pathname)
   %READBESTSNOWDEPTH Site-aware snow-depth selector for ESM-SnowMIP obs.
   %
   %  [values, source] = icemodel.verification.setup.readBestSnowDepth(pathname)
   %
   %  ESM-SnowMIP boreal forest sites (oas, obs, ojp) report snd_gap_auto /
   %  snd_gap1_auto in the canopy gap rather than snd_auto. This selector
   %  prefers, in order:
   %    snd_auto > snd_gap_auto > snd_gap1_auto > snd_man
   %  so the resulting series is the most representative open-area
   %  snow-depth time series available for the site.
   %
   %  Returns
   %    values : double column [m]
   %    source : string identifying which channel was used
   %
   % See also: icemodel.verification.setup.readObsChannel

   info = ncinfo(pathname);
   names = string({info.Variables.Name});
   candidates = ["snd_auto", "snd_gap_auto", "snd_gap1_auto", "snd_man"];
   for c = candidates
      if any(names == c)
         values = icemodel.verification.setup.readObsChannel(pathname, c);
         source = c;
         return
      end
   end
   error('icemodel:verification:setup:readBestSnowDepth:noChannel', ...
      'no usable snow-depth channel in %s', pathname);
end
