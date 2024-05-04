function mesh(Z, dz, growth_factor)
   % Plot the nodes and edges of an exponential grid.
   %
   % Parameters:
   %   Z (float): Total depth of the grid.
   %   dz (float): Initial control volume thickness.
   %   growth_factor (float): Factor by which CV thickness increases each layer.
   %
   % Example usage
   % icemodel.plot.mesh(10, 1, 1.1);

   arguments
      Z (1, 1) {mustBeNumeric} = 20
      dz (1, 1) {mustBeNumeric} = 0.04
      growth_factor (1, 1) {mustBeNumeric} = 1.1
   end

   % Initialize arrays and variables
   z = 0;         % Start at the top
   n = 1;         % Layer counter
   Z_tot = 0;     % Total thickness covered
   nodes = [];    % Initialize node array
   edges = 0;     % Start with the top edge

   % Compute positions of edges and nodes
   while Z_tot < Z
      current_dz = dz * growth_factor^(n-1);
      Z_tot = Z_tot + current_dz;

      if Z_tot > Z
         current_dz = current_dz - (Z_tot - Z);
         Z_tot = Z;
      end

      z = z + current_dz;
      nodes(n) = z - current_dz / 2; %#ok<*AGROW>
      edges(n + 1) = z;

      n = n + 1;
   end

   % Plot nodes as filled circles
   figure; hold on
   plot(ones(size(nodes)), nodes, 'bo', 'MarkerFaceColor', 'b', ...
      'MarkerEdgeColor', 'none'); % Interior
   plot(1, 0, 'ro', 'MarkerFaceColor', 'r', ...
      'MarkerEdgeColor', 'none'); % Top Boundary
   plot(1, Z, 'ro', 'MarkerFaceColor', 'r', ...
      'MarkerEdgeColor', 'none'); % Bottom Boundary

   % Plot edges as horizontal lines
   for n = 1:length(edges)
      line([0.8, 1.2], [edges(n), edges(n)], 'Color', 'k');
   end

   % Formatting
   title('Exponential Grid');
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
   % figure; hold on
   % plot(x, z_node_bc, 'o', ...
   %    'MarkerEdgeColor', 'none');
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
   % plot(ones(size(z_spect)), z_spect, 'o', ...
   %    'MarkerEdgeColor', 'none');
   % for n = 1:numel(z_walls)
   %    plot([x(n)/2 x(n)+x(n)/2], [z_walls(n) z_walls(n)], 'k', 'LineWidth', 1)
   % end
   % set(gca, 'YDir', 'reverse')
   %
   %
   % plot([ones(size(z_walls)) ones(size(z_walls))]', [z_walls z_walls]')
end
