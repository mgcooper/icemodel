function assertNewReleaseBaselineTarget(output_file, overwrite)
   %ASSERTNEWRELEASEBASELINETARGET Reject an existing immutable release file.

   arguments
      output_file (1, 1) string
      overwrite (1, 1) logical = false
   end

   if ~isfile(output_file)
      return
   end

   if overwrite
      overwrite_note = " overwrite=true cannot replace an immutable release.";
   else
      overwrite_note = "";
   end
   error('icemodel:test:releaseBaselineImmutable', ...
      'Release baseline already exists: %s.%s', output_file, overwrite_note)
end
