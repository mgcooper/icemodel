function plotncfile(f, gridcells)

   info = ncparse(f);
   data = ncreaddata(f);

   if ismember("Tice", info.Name)

      plotice2(data, info, gridcells);

   elseif ismember("Tsfc", info.Name)

      plotice1(data, info, gridcells);
   end
end

function plotice2(data, info, gridcells)

   depth = data.depth;
   tice = data.Tice;

   % Unless/until ncreaddata orientation of data changes, the data will be:
   % [gridcell x time x depth], thus below we take the first dimension to get
   % one grid cell and average down the time dimension to get a depth profile.

   % Extract the grid cells then take the temporal average
   tice = tice(gridcells, :, :);

   tice_avg = squeeze(mean(tice, 1));

   % Plot
   figure
   plot(tice_avg, depth)
   set(gca, 'YDir', 'Reverse')
   xlabel("T [oC]")
   ylabel('depth [m]')
   % title('annual average')

   % Plot each grid cell if fewer than 15
   if numel(gridcell) < 15
      % or mabye call this from a loop

   end


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
