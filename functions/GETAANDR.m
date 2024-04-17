function [a, r] = GETAANDR(bulkcoefs, albedo)
   %GETAANDR Get the bulk absorption and reflection coefficients at each layer
   a = (1.0 - albedo) ./ (1.0 + albedo) .* bulkcoefs;
   r = 2.0 .* albedo .* bulkcoefs ./ (1.0 - albedo ^ 2);
end
