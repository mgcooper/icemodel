function merged = preferPrimary(primary, fallback)
   %PREFERPRIMARY Coalesce two same-shape series, keeping primary where finite.
   %
   %  merged = icemodel.verification.setup.preferPrimary(primary, fallback)
   %
   %  Returns the primary series, with non-finite entries replaced by the
   %  corresponding fallback entries when those are finite. Used by the
   %  ESM-SnowMIP builders to merge automatic / manual gauge channels
   %  (e.g. snd_auto with snd_man fallback).

   merged = primary;
   replace = ~isfinite(merged) & isfinite(fallback);
   merged(replace) = fallback(replace);
end
