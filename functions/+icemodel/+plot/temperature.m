function temperature(T, dz)
   % need this in a function, and an anonymous one to create z_therm:

   znodes = @(dz) cumsum(dz) - dz / 2;

   figure;
   plot(T, znodes(dz)); set(gca, 'YDir', 'Reverse')
end
