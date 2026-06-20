function case_types = casetype()
   %CASETYPE Return the supported snow-verification case types.
   %
   %  case_types = icemodel.verification.namelists.casetype()
   %
   % Outputs
   %  case_types   Supported manifest case_type values. These classify the
   %               verification case itself, not the source dataset family.
   %
   % Role
   %  Canonical manifest classification list. Case types answer "what kind of
   %  verification case is this?", for example field site or synthetic process.

   % Keep this explicit so adding a new group is a deliberate schema change.
   % "firn_observational" is the soft-gated firn evaluation case type (reports
   % diagnostic metrics, no hard PASS/FAIL tolerance); "firn_analytical" is
   % deferred with the Meyer-Hewitt namespace and is intentionally absent.
   case_types = ["esm_site"; "synthetic_process"; "firn_observational"];
end
