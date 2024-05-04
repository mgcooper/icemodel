function dwavl = GETDWAVL(wavelength, nwavl)
   %GETDWAVL Compute first derivative of wavelength
   %
   %#codegen

   dwavl = zeros(1, nwavl);
   dwavl(1) = 2.0 * (wavelength(1, 2) - wavelength(1, 1));
   for k = 2:nwavl-1
      dwavl(k) = (wavelength(1, k+1) - wavelength(1, k-1)) / 2.0;
   end
   dwavl(nwavl) = 2.0 * (wavelength(1, nwavl) - wavelength(1, nwavl-1));
   % [um] (micrometers)
end
