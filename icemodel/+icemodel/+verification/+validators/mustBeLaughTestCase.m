function mustBeLaughTestCase(case_ids)
   %MUSTBELAUGHTESTCASE Validate cases against the Laugh-Tests namelist.
   %
   %  icemodel.verification.validators.mustBeLaughTestCase(case_ids)
   %
   % MATLAB arguments blocks cannot call a dynamic namelist inside mustBeMember.
   % This validator keeps that lookup canonical while remaining arguments-safe.

   valid = icemodel.verification.namelists.caseid("laugh_tests");
   bad = setdiff(reshape(string(case_ids), 1, []), valid, 'stable');
   if ~isempty(bad)
      error('icemodel:verification:validators:mustBeLaughTestCase', ...
         'unknown Laugh-Tests case %s. Valid: %s', ...
         bad(1), strjoin(valid, ', '))
   end
end
