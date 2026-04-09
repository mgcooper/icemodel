function Qc = conductive_heat_flux(k_eff, T, dz, T_sfc, level)
   %CONDUCTIVE_HEAT_FLUX Compute conductive heat flux, positive into the surface.
   %
   %  Qc = icemodel.surface.conductive_heat_flux(k_eff, T, dz, T_sfc)
   %  Qc = icemodel.surface.conductive_heat_flux(k_eff, T, dz, T_sfc, level)
   %
   % level = 1 (default): conduction across the top half-cell boundary
   %   Qc = k_eff(1) * (T(1) - T_sfc) / (dz(1) / 2)
   %
   % level = 2: conduction across the interior cell boundary between layers 1 and 2
   %   Qc = (k_eff(1) + k_eff(2)) / 2 * (T(2) - T(1)) / ((dz(1) + dz(2)) / 2)
   %
   % Units: [W m-2] = [W m-1 K-1] * [K] / [m]
   %
   % See also: icemodel.surface.solve_surface_energy_balance,
   %           icemodel.surface.diagnose_melt_freeze_energy
   %
   %#codegen

   if nargin < 5
      level = 1;
   end

   if level == 1
      % Conduction from layer 1 into the surface (across a 1/2 cv)
      Qc = k_eff(1) * (T(1) - T_sfc) / (dz(1) / 2);
   elseif level == 2
      % Conduction from layer 2 into layer 1 (across a full cv)
      Qc = (k_eff(1) + k_eff(2)) / 2.0 * (T(2) - T(1)) / (dz(1) + dz(2)) / 2.0;
   end
end
