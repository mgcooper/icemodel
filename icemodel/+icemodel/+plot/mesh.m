function plotmesh(Z, dz, growth_factor)
   % Plot the nodes and edges of an exponential grid.
   %
   % Parameters:
   %   Z (float): Total depth of the grid.
   %   dz (float): Initial control volume thickness.
   %   growth_factor (float): Factor by which CV thickness increases each layer.

   % Example usage
   % plot_exponential_grid(10, 1, 1.1);


   % Initialize arrays and variables
   z = 0;         % Start at the top
   n = 1;         % Layer counter
   total_thickness = 0; % Total thickness covered
   nodes = [];    % Initialize node array
   edges = [0];   % Start with the top edge

   % Compute positions of edges and nodes
   while total_thickness < Z
      current_dz = dz * growth_factor^(n-1);
      total_thickness = total_thickness + current_dz;

      if total_thickness > Z
         current_dz = current_dz - (total_thickness - Z);
         total_thickness = Z;
      end

      z = z + current_dz;
      nodes(n) = z - current_dz / 2;
      edges(n + 1) = z;

      n = n + 1;
   end

   % Plotting
   figure;

   % Plot nodes as filled circles
   plot(ones(size(nodes)), nodes, 'bo', 'MarkerFaceColor', 'b'); % Interior Nodes
   hold on;
   plot(1, 0, 'ro', 'MarkerFaceColor', 'r'); % Top Boundary Node
   plot(1, Z, 'ro', 'MarkerFaceColor', 'r'); % Bottom Boundary Node

   % Plot edges as horizontal lines
   for i = 1:length(edges)
      line([0.8, 1.2], [edges(i), edges(i)], 'Color', 'k');
   end

   % Formatting
   title('Exponential Grid Visualization');
   xlabel('Grid');
   ylabel('Depth');
   set(gca, 'YDir', 'reverse'); % Invert y-axis
   xlim([0.5, 1.5]);
   ylim([0, Z]);
   legend('Interior Nodes', 'Boundary Nodes', 'Location', 'best');
   grid on;
   hold off;


   % %% Use this if within CVMESH
   %
   % x = ones(size(z_node_bc));
   %
   % figure;
   % plot(x, z_node_bc, 'o'); hold on
   % formatPlotMarkers('markersize', 6)
   %
   % for n = 1:numel(z_edge)
   %    plot([x(n)/2 x(n)+x(n)/2], [z_edge(n) z_edge(n)], 'k', 'LineWidth', 1)
   % end
   %
   % set(gca,'YDir', 'reverse', 'YScale', 'linear')
   % % set(gca,'YDir', 'reverse', 'YScale', 'log')
   %
   %
   % %% Use this if within UPDATEEXTCOEFS
   %
   % x = ones(size(z_walls));
   %
   % figure
   % hold on
   % plot(ones(size(z_spect)), z_spect, 'o')
   % for n = 1:numel(z_walls)
   %    plot([x(n)/2 x(n)+x(n)/2], [z_walls(n) z_walls(n)], 'k', 'LineWidth', 1)
   % end
   %
   % formatPlotMarkers('markersize', 6)
   % set(gca, 'YDir', 'reverse')
   %
   %
   % plot([ones(size(z_walls)) ones(size(z_walls))]', [z_walls z_walls]')

end

