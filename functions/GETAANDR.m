%--------------------------------------------------------------------------
%   get the bulk absorption and reflection coefficients at each layer
%--------------------------------------------------------------------------

function [a,r] = GETAANDR(bulkcoefs,albedo)
%--------------------------------------------------------------------------
    a = (1.0 - albedo) ./ (1.0 + albedo) .* bulkcoefs;
    r = 2.0 .* albedo .* bulkcoefs ./ (1.0 - albedo^2);