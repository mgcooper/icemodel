function bulkextcoefs(bulkcoefs)
   %icemodel.plot.bulk_extinction_coefficients
   %
   % Call this from within icemodel.radiation.bulk_extinction_coefficients

   figure;
   try
      semilogx(bulkcoefs, z_spect);
   catch
      semilogx(bulkcoefs(1:numel(z_spect)), z_spect);
   end
   set(gca,'YDir', 'reverse', 'XDir', 'reverse');
   xlabel('\eta (z)')
   ylabel('z')
   % copygraphics(gcf)
end
