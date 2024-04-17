function varargout = plotice2(ice2, varname, iter, zdepth, dz, dt)
   %PLOTICE2 Plot 2-d icemodel data

   if isnumeric(varname)
      varname = inputname(1);
   end

   if nargin < 4
      zdepth = 2;
   end
   if nargin < 5
      try
         dz = ice2.Z(2) - ice2.Z(1);
      catch
         dz = 0.04;
      end
   end
   if nargin < 6
      dt = 3600;
   end

   zidx = ceil(zdepth / dz(1));
   vmat = ice2.(varname)(1:zidx, 1:iter);
   tmat = dt * (1:iter);
   zmat = 0:dz(1):zdepth;

   if numel(zmat) ~= size(vmat, 1)
      zmat = 0:dz(1):zdepth-dz(1);
   end

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

   if nargout == 1
      varargout{1} = H;
   end
end
