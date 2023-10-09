function [solar,total_solar] = GETSOLAR(solar,nvalues,wavelength,dwavelen)
%GETSOLAR Interpolate the spectral solar radiation to the 118 bands of mie.dat,
%and then integrate it to compute total_solar.
isolarvals = size(solar,1);
wavel_tmp = solar(:,1);
solar_tmp = solar(:,2);

% Generate a dummy downward solar spectrum using the above _tmp
%   data strings, interpolating to the wavelengths of interest.
temp = nan(nvalues,1);
for k=1:nvalues
   x = wavelength(1,k);
   for i=1:isolarvals-1
      if (x>wavel_tmp(i))
         icount = i;
      end
   end
   x1 = wavel_tmp(icount);
   x2 = wavel_tmp(icount+1);
   y1 = solar_tmp(icount);
   y2 = solar_tmp(icount+1);

   temp(k) = y1 + (x - x1) * (y2 - y1)/(x2 - x1);
end

% mgc this is different than glen's code
solar = temp';

% Integrate the solar radiation.
total_solar = 0.0;
for k=1:nvalues
   total_solar = total_solar + solar(k) * dwavelen(k);
end
