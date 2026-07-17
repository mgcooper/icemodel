function text = formatManifestTime(value)
   %FORMATMANIFESTTIME Serialize manifest timestamps with explicit clock time.
   %
   %  text = icemodel.verification.setup.formatManifestTime(value)
   %
   % Manifest periods are read by humans and by staging helpers that compare
   % windows across preserved artifacts. MATLAB's default string conversion can
   % omit the clock for midnight values, so this helper writes every finite
   % bound as yyyy-MM-dd HH:mm:ss while preserving blank / NaT values for
   % intentionally unbounded periods.

   if isstring(value) || ischar(value)
      if all(strlength(string(value)) == 0)
         text = '';
         return
      end
   elseif isdatetime(value) && all(isnat(value))
      text = '';
      return
   end

   % Normalize to the repository's UTC convention before assigning a display
   % format. The returned char is format-stable regardless of MATLAB defaults.
   value = icemodel.verification.setup.ensureUtc(value);
   value.Format = 'yyyy-MM-dd HH:mm:ss';
   text = char(string(value));
end
