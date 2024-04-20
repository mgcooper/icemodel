% need an accounting of all places where indexing is used
%
% defdatavars
% defdimvars
% getvardata
% writedims - ok b/c it queries the vid first
% writeice1 - ok b/c it queries the vid first
%
% need to see indexing like this:
% fields = string(fieldnames(ncatts));
% for thisfield = fields(:)'
%   netcdf.putAtt(ncid, varid, thisfield, ncatts.(thisfield));
% end
%
% not this:
% for v = 1:numel(vars)
%   thisvar = vars{v};
%   thisvid = netcdf.defVar(ncid, thisvar, xtype, dimid);
% end
%
%
% if Z, dz are default 0, then option "getDimsFromData" should work for ice2
% the reason it doesn't is b/c