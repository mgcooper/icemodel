function V = interpRcm(x, y, v, xq, yq, kwargs)
   %INTERPRCM Interpolate RCM grid-cell data to query points.
   %
   %  V = icemodel.forcing.helpers.interpRcm(x, y, v, xq, yq)
   %  V = ... interpRcm(_, method="natural")
   %
   % Scattered interpolation of regional-climate-model (RCM) grid-cell
   % data onto query points (a site location, a set of points, or the
   % interpolation grid covering a catchment polygon). The
   % triangulation is built once and reused across timesteps.
   %
   % Inputs
   %  x, y   - projected coordinates of the RCM grid cells (n cells)
   %  xq, yq - query point coordinates in the same projection (m points)
   %  v      - data, cells x time (time x cells is auto-oriented when
   %           unambiguous)
   %  method - (name-value) "natural" (default), "linear", or "nearest"
   %
   % Outputs
   %  V - interpolated data, time x m (timeseries down the columns)
   %
   % See also: scatteredInterpolant, icemodel.forcing.helpers.psnProjection

   arguments
      x (:, 1) double
      y (:, 1) double
      v (:, :) double
      xq (:, 1) double
      yq (:, 1) double
      kwargs.method (1, 1) string {mustBeMember(kwargs.method, ...
         ["natural", "linear", "nearest"])} = "natural"
   end

   assert(numel(x) == numel(y), 'X and Y must have equal lengths')
   assert(numel(xq) == numel(yq), 'XQ and YQ must have equal lengths')

   % Orient v as cells down the rows, time across the columns.
   if size(v, 1) ~= numel(x)
      v = v.';
   end
   assert(size(v, 1) == numel(x), ...
      'V must be cells x time, with cells matching X and Y')

   % Interpolate each timestep, reusing the triangulation by updating
   % the interpolant values (much faster than rebuilding per step).
   n_time = size(v, 2);
   V = nan(numel(xq), n_time);
   F = scatteredInterpolant(x, y, v(:, 1), kwargs.method);
   V(:, 1) = F(xq, yq);
   for n = 2:n_time
      F.Values = v(:, n);
      V(:, n) = F(xq, yq);
   end

   % Return timeseries down the columns.
   V = V.';
end
