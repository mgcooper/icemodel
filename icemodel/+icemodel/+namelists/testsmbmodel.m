function list = testsmbmodel()
%TESTSMBMODEL Return the supported formal test smbmodel selectors.
%
%  list = icemodel.namelists.testsmbmodel()

   list = ["all"; icemodel.namelists.smbmodel()];
end
