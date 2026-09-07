function [d_pevp, pevp] = potential_surface_vapor_demand(Qe, dt, dz)
   %POTENTIAL_SURFACE_VAPOR_DEMAND Convert latent heat to surface demand.
   %
   %  [d_pevp, pevp] = ...
   %     icemodel.kernels.potential_surface_vapor_demand(Qe, dt, dz)
   %
   % QE is the latent heat flux [W m-2]. PEVP is its liquid-water-
   % equivalent velocity [m s-1], and D_PEVP is the corresponding top-cell
   % liquid-water-equivalent fraction over DT [s] and DZ [m]. Positive is
   % condensation or deposition; negative is evaporation or sublimation.
   % The surface exchange partitions this energy demand between liquid at
   % Lv and ice at Ls.
   %
   % Definition:
   %
   %   pevp   = Qe / (Lv * ro_liq)
   %   d_pevp = pevp * dt / dz
   %          = Qe / (Lv * ro_liq) * dt / dz
   %
   % The first conversion changes latent heat flux [J s-1 m-2] to
   % liquid-water-equivalent velocity [m s-1]. Time integration and division
   % by top-cell thickness produce the dimensionless fraction demand:
   %
   %   [m s-1] = [J s-1 m-2] / ([J kg-1] * [kg m-3])
   %   [-] = [m s-1] * [s] / [m]
   %
   % The ice-fraction equivalent used by sublimation accounting is
   %
   %   d_psbl = d_pevp * (Lv * ro_liq) / (Ls * ro_ice).
   %
   % See also: icemodel.surface.potential_surface_vapor_demand,
   %  icemodel.surface.potential_surface_vapor_exchange
   %
   %#codegen

   persistent Lv ro_liq
   if isempty(Lv)
      [Lv, ro_liq] = icemodel.physicalConstant("Lv", "ro_liq");
   end

   pevp = Qe / (Lv * ro_liq);
   d_pevp = pevp * dt / dz;
end
