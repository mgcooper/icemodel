function list = rcmMetSources()
   %RCMMETSOURCES Verification RCM labels that currently write met files.
   %
   %  list = icemodel.verification.namelists.rcmMetSources()
   %
   % This is intentionally narrower than rcmsources: RACMO 2.3p3 is staged as
   % userdata/Data only because the available subsurface product does not carry
   % the near-surface met state channels required by validatemet.

   list = ["mar", "merra"];
end
