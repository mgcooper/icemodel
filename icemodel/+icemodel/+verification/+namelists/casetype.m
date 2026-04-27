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
   case_types = ["esm_site"; "synthetic_process"];
end
