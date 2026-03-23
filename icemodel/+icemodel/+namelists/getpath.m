function list = getpath()
%GETPATH Return the supported icemodel.getpath path kinds.
%
%  list = icemodel.namelists.getpath()

   list = ["data"; "demo"; "input"; "eval"; "userdata"; ...
      "output"; "restart"; "test"];
end
