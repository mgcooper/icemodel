function case_ids = laughtests()
   %LAUGHTESTS Canonical Laugh-Tests case-id namelist.
   %
   %  case_ids = icemodel.verification.namelists.laughtests()
   %
   % Outputs
   %  case_ids   String column of runnable case ids in the "laugh_tests"
   %             dataset family.
   %
   % Role
   %  Mirrors icemodel.verification.namelists.snowmipsite for the
   %  Laugh-Tests dataset family so case-id resolution stays uniform
   %  across families. Currently only colbeck1976 is staged; new
   %  Laugh-Tests cases are added by extending this list.

   case_ids = "colbeck1976";
end
