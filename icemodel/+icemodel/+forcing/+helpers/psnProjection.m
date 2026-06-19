function proj = psnProjection()
   %PSNPROJECTION Polar stereographic north projection used by the builders.
   %
   %  proj = icemodel.forcing.helpers.psnProjection()
   %
   % Returns the WGS 84 / NSIDC Sea Ice Polar Stereographic North
   % projection (EPSG:3413, latitude of true scale 70N, central meridian
   % 45W) as a projcrs object. This is the projection the legacy runoff
   % builders loaded from projsipsn.mat. Use projfwd / projinv with the
   % returned object to convert geographic lat/lon to projected x/y.
   %
   % See also: projcrs, projfwd, projinv

   proj = projcrs(3413);
end
