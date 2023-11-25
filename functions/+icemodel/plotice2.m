function H = plotice2(ice2, varname, iter, zdepth, dz, dt)
   %PLOTICE2 Plot 2-d icemodel data 

   if isnumeric(varname)
      varname = inputname(1);
   end

   zidx = ceil(zdepth / dz(1));
   vmat = ice2.(varname)(1:zidx, 1:iter);
   tmat = dt * (1:iter);
   zmat = 0:dz:zdepth-dz(1);

   maxfig;
   H = pcolor(tmat / 86400, zmat, vmat);
   shading flat
   
   c = colorbar(gca, 'eastoutside');
   setcolorbar(c, 'Label', varname)
   xlabel({'time (days)', ''})
   ylabel('depth (m)')
   
   set(gca, ...
      'YDir', 'reverse', ...
      'PositionConstraint', 'outerposition', ...
      'TickLength', [0.01 0.01]);
   
   % From earlier versions
   % H = pcolor(tmat / 86400, zmat, vmat);
   % H.EdgeColor = 'none';
   % 
   % setcolorbar(c, 'Title', varunits)
   % title(varname);
   
end
