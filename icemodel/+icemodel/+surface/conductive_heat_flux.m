function [Qc, dQc_dT_sfc] = conductive_heat_flux(k_eff, T, dz, T_sfc)
   %CONDUCTIVE_HEAT_FLUX Conductive heat flux into the surface and derivative.
   %
   %  Qc = icemodel.surface.conductive_heat_flux(k_eff, T, dz, T_sfc)
   %  [Qc, dQc_dT_sfc] = icemodel.surface.conductive_heat_flux( ...
   %     k_eff, T, dz, T_sfc)
   %
   % Computes the conductive heat flux from ice layer 1 into the surface
   % across the top half-control-volume boundary:
   %
   %   Qc = k_eff(1) * (T(1) - T_sfc) / (dz(1) / 2)   [W m^-2]
   %
   % The partial derivative with respect to T_sfc is constant and negative:
   %
   %   dQc/dT_sfc = -k_eff(1) / (dz(1) / 2)            [W m^-2 K^-1]
   %
   % solve_surface_temperature includes the derivative in the Newton-Raphson
   % Jacobian for the Dirichlet solve. For the Robin solve, conduction enters
   % through the top-node finite-difference equation in
   % icemodel.column.assemble_enthalpy_system rather than this derivative.
   %
   % See also: icemodel.surface.solve_surface_temperature,
   %           icemodel.surface.diagnose_melt_freeze_energy
   %
   %#codegen

   % Qc equals the enthapy solve's top face term a1 (see
   % icemodel.column.assemble_enthalpy_system). Both use vapor-free k_eff(1)
   % because face 1 is closed for vapor. Keep the SEB and column solve
   % consistent by passing the solver's vapor-free k_eff to this function.
   Qc = k_eff(1) * (T(1) - T_sfc) / (dz(1) / 2);

   if nargout > 1
      dQc_dT_sfc = -k_eff(1) / (dz(1) / 2);
   end
end
