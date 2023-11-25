function plotcvmix(f_ice1, f_liq1, f_ice2, f_liq2, f_air1, f_air2, dz)
   %PLOTCVMIX Plot control volume mixture 

   if nargin < 5
      f_air1 = 1 - f_liq1 - f_ice1;
      f_air2 = 1 - f_liq2 - f_ice2;
   end
   
   if nargin < 7
      dz = 1.0;
   end

   X = [1 2 3];
   Y = [ 
      f_ice1 f_liq1 f_air1;               % cv 1
      f_ice2 f_liq2 f_air2;               % cv 2
      f_ice1+f_ice2 f_liq1+f_liq2 0;      % liq + ice, no air
      ];

   Y = dz * Y;
   
   % Make the figure
   figure;
   bar(X, Y,'stacked'); hold on;
   legend('ice','liq','air');
   ylabel('fraction (-)');
end
