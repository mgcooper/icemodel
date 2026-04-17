function [I0, solar, solar_dwavel] = get_solar_spectrum(solar_spectrum, wavel, dwavel)
   %get_solar_spectrum Interpolate the reference solar spectrum to the model grid.
   %
   %  [I0, solar, solar_dwavel] = icemodel.radiation.get_solar_spectrum( ...
   %     solar_spectrum, wavel, dwavel)
   %
   % solar_spectrum is a two-column array, the first column is wavelength, the
   % second column is solar irradiance.
   %
   % The spectral model works with one fixed prototype solar spectrum. I0 is the
   % corresponding integrated incoming shortwave over that prototype spectrum.
   %
   % The model ships optical properties (mie.dat) on a defined 118-band
   % wavelength grid. The solar spectrum is interpolated to these bands here.
   %
   %#codegen

   % Interpolate the tabulated solar spectrum onto the model wavelength grid.
   solar = interp1(solar_spectrum(:, 1), solar_spectrum(:, 2), ...
      wavel, 'linear');
   solar = solar(:).';

   % Integrate the interpolated spectrum with the configured quadrature
   % weights to recover the total incoming prototype flux.
   solar_dwavel = solar .* dwavel;
   I0 = sum(solar_dwavel);
end
