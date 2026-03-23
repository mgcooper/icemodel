function list = smbmodel(group)
   %SMBMODEL Return supported smbmodel names by group.
   %
   %  list = icemodel.namelists.smbmodel()
   %  list = icemodel.namelists.smbmodel("test")

   arguments
      group (1, :) string {mustBeMember(group, ["all", "test"])} = "all"
   end

   % Return the full public smbmodel list by default.
   switch group
      case "all"
         list = ["icemodel"; "skinmodel"; "firnmodel"];

      case "test"
         list = ["icemodel"; "skinmodel"];
   end
end
