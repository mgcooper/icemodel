function list = physicalconstant()
%PHYSICALCONSTANT Return the supported physical constant names.
%
%  list = icemodel.namelists.physicalconstant()

   constants = icemodel.physicalConstant();
   list = string(fieldnames(constants));
end
