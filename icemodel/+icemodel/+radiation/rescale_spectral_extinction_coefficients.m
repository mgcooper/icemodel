function k_ext = rescale_spectral_extinction_coefficients(k_ext, kabs, kice, wavel)
   %rescale_spectral_extinction_coefficients Scale extinction coefficients by an absorption profile.
   %
   %  k_ext = icemodel.radiation.rescale_spectral_extinction_coefficients( ...
   %     k_ext, kabs, kice, wavel)
   %  Scales the theoretical spectral extinction coefficients by a
   %  user-provided absorption coefficient profile from 0.3-0.9 um.
   %
   % Use Eq. 15 of Warren et al. (2006) to scale the asymptotic two-stream
   % spectral extinction coefficients by a user-defined absorption spectrum for
   % ice or snow.
   %
   % Stephen G. Warren, Richard E. Brandt, and Thomas C. Grenfell
   % "Visible and near-ultraviolet absorption spectrum of ice from
   % transmission of solar radiation into snow," Applied Optics,
   % Vol. 45, No. 21.
   %
   %#codegen

   % The spectral bands used here do not include 0.60 um, so interpolate the
   % reference values between the neighboring bands at indices 14 and 15. This
   % could be hard-coded from Warren et al. (2008), but keeping the
   % interpolation here preserves the direct relationship to the active grid.
   dwavel = wavel(15) - wavel(14);
   dkabs = kice(15) - kice(14);
   wavel0 = 0.60;
   kice0 = kabs(14) + (dkabs / dwavel) * (wavel0 - wavel(14));

   % Compute the matching reference extinction coefficient at 0.60 um.
   dkext = k_ext(15) - k_ext(14);
   kext0 = k_ext(14) + (dkext / dwavel) * (wavel0 - wavel(14));

   % Apply Eq. 15:
   kext = kext0 * sqrt(kabs ./ kice0);

   % Merge the scaled values from 0.3-0.9 um with the rest of the spectrum.
   % k_ext_scaled = [kext(1:24) k_ext(25:end)];
   k_ext = [kext(1:24), k_ext(25:end)];

   % Use this to see the result.
   % figure;
   % plot(wavel, k_ext); hold on; set(gca, 'YScale', 'log');
   % plot(wavel, k_ext_scaled);
   % legend('reference', 'user defined');
end
