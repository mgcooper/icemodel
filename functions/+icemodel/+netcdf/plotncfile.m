function plotncfile(f)

   info = ncparse(f);
   data = ncreaddata(f);

   if ismember("Tice", info.Name)

      plotice2(data, info);

   elseif ismember("Tsfc", info.Name)

      plotice1(data, info);
   end
end

function plotice2(data, info)
   
   depth = data.depth;
   tice = data.Tice;
   
   % Unless/until ncreaddata orientation of data changes, the data will be:
   % [gridcell x time x depth], thus below we take the first dimension to get
   % one grid cell and average down the time dimension to get a depth profile.
   
   % Extract one grid cell then take the temporal average
   tice = mean(squeeze(tice(1, :, :)), 1);
   
   % Plot 
   figure
   plot(tice, depth)
   set(gca, 'YDir', 'Reverse')
   xlabel("T [oC]")
   ylabel('depth [m]')
   
   
   % might be able to levereage this
   % icemodel.plotice2(ice2_38, 'Tice', 8760, 20, 0.04)
   
end


function plotice1(data, info)

   x = data.x_easting;
   y = data.y_northing;

   % Plot all 1d vars
   vars = {'Tsfc', 'melt', 'freeze', 'runoff', 'cond', 'subl'};
   figontop
   for n = 1:numel(vars)
      v = data.(vars{n});
      scatter(x, y, 20, mean(v, 2), 'filled')
      colorbar
      title(vars{n})
      pause
      clf
   end
end
