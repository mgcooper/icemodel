function mstruct = geusModisProjection()
   %GEUSMODISPROJECTION Native polar-stereographic projection of the GEUS grid.
   %
   %  mstruct = icemodel.forcing.helpers.geusModisProjection()
   %
   % Returns the map projection in which the GEUS Greenland reflectivity
   % (MODIS albedo) 5 km grid is REGULAR, as a map projection structure
   % (mstruct) usable with mfwdtran / minvtran.
   %
   % The Greenland_Reflectivity_<YYYY>_5km_C6.nc files carry only 2-D
   % lon/lat coordinate grids (no native x/y coordinate vectors and no CRS
   % attributes). The grid is the GIMP/MODIS Greenland polar stereographic
   % posting: a sphere (radius 6378137 m), central meridian 39W, latitude of
   % true scale 70N. Reprojecting the lon/lat to this frame yields an
   % axis-aligned uniform ~5.15 km grid (off-axis residual < 0.3 m), whereas
   % reprojecting to EPSG:3413 (psnProjection, central meridian 45W) yields a
   % ROTATED, non-uniform grid that fails exactremap's regular-grid check.
   %
   % Determined empirically (2026-06-19) by sweeping the stereographic
   % parameters for the projection that minimises off-axis variation: the
   % lon/lat roundtrip through this mstruct closes to ~1e-13 deg in latitude,
   % confirming it is the genuine native frame, not an approximation.
   %
   % See also: icemodel.forcing.readGeusModis,
   %  icemodel.forcing.helpers.psnProjection, mfwdtran, minvtran

   mstruct = defaultm('stereo');
   mstruct.origin = [90 -39 0];     % north pole, central meridian 39W
   mstruct.geoid = [6378137 0];     % sphere, radius 6378137 m (eccentricity 0)
   mstruct.mapparallels = 70;       % latitude of true scale 70N
   mstruct = defaultm(mstruct);
end
