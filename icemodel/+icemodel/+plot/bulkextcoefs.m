function bulkextcoefs(bulkcoefs)

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
