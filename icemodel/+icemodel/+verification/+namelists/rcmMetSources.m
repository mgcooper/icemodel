function list = rcmMetSources()
   %RCMMETSOURCES Verification RCM labels that currently write met files.
   %
   %  list = icemodel.verification.namelists.rcmMetSources()
   %
   % This list is deliberately narrower than rcmsources. RACMO 2.3p3 stages as
   % userdata/Data only, because the available subsurface product does not
   % carry the near-surface met state channels that validatemet requires.

   list = ["mar", "merra"];
end
