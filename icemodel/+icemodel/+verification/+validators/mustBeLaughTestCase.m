function mustBeLaughTestCase(case_id)
   %MUSTBELAUGHTESTCASE Validate one case against the Laugh-Tests namelist.
   %
   %  icemodel.verification.validators.mustBeLaughTestCase(case_id)
   %
   % MATLAB arguments blocks cannot call a dynamic namelist inside mustBeMember.
   % This validator keeps that lookup canonical while remaining arguments-safe.

   valid = icemodel.verification.namelists.caseid("laugh_tests");
   if ~ismember(case_id, valid)
      error('icemodel:verification:validators:mustBeLaughTestCase', ...
         'unknown Laugh-Tests case %s. Valid: %s', ...
         case_id, strjoin(valid, ', '))
   end
end
