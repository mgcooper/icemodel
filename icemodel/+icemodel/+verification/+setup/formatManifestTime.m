function text = formatManifestTime(value)
   %FORMATMANIFESTTIME Serialize manifest timestamps with explicit clock time.
   %
   %  text = icemodel.verification.setup.formatManifestTime(value)
   %
   % People read manifest periods, and so do staging helpers that compare
   % windows across preserved artifacts. MATLAB's default string conversion
   % can omit the clock for midnight values. This helper therefore writes
   % every finite bound as yyyy-MM-dd HH:mm:ss. It keeps blank and NaT values
   % for periods that are unbounded on purpose.

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
