function [Sc, chi] = shortwave_source_term(Qsi, albedo, I0, dz_spect, tau_N, ...
      tau_S, solar_dwavel, dz_therm, ro_sno, z_nodes_therm, z_nodes_spect, ...
      z_edges_spect, k_bulk_lookup)
   %shortwave_source_term Solve the spectral shortwave source term.
   %
   % [Sc, chi] = icemodel.column.shortwave_source_term(Qsi, albedo, I0, ...
   %    dz_spect, tau_N, tau_S, solar_dwavel, dz_therm, ro_sno, z_nodes_therm, ...
   %    z_nodes_spect, z_edges_spect, k_bulk_lookup)
   %
   % The spectral model solves for the net shortwave flux profile on the
   % spectral control-volume grid, then collapses that absorbed-flux profile to
   % the thermal grid used by the enthalpy solver.
   %
   % If the K_BULK_LOOKUP argument is an empty struct, it switches the
   % bulk-extinction step from the density lookup approximation to the exact
   % solution. All other parts of the spectral solve remain shared.
   %
   %#codegen

   % Dark/no-sun early exit.
   if Qsi < 1e-3
      Sc = 0.0 * ro_sno;
      chi = 1.0;
      return
   end

   % Map the thermal-grid density profile onto the spectral grid.
   ro_sno_spect = interp1(z_nodes_therm, ro_sno, z_nodes_spect, ...
      'nearest', 'extrap');

   % Build the bulk extinction coefficients with either the exact transform or
   % the lookup-table shortcut.
   if isempty(k_bulk_lookup)
      k_bulk = icemodel.radiation.bulk_extinction_coefficients( ...
         dz_spect, ro_sno_spect, tau_N, tau_S, solar_dwavel);
   else
      k_bulk = icemodel.radiation.bulk_extinction_coefficients_lookup( ...
         ro_sno_spect, k_bulk_lookup);
   end

   % Solve the two-stream system for the net flux at each interface.
   Qnet = icemodel.radiation.solvetwostream(I0, albedo, k_bulk, z_edges_spect);

   % Collapse the spectral-grid net flux to the thermal-grid source term and
   % compute the chi partition used by the SEB coupling.
   [Sc, chi] = spectralNetFluxToSourceTerm(Qsi, I0, albedo, Qnet, ...
      dz_spect, dz_therm);
end

function [Sc, chi] = spectralNetFluxToSourceTerm(Qsi, I0, albedo, Qnet, ...
      dz_spect, dz_therm)
   %SPECTRALNETFLUXTOSOURCETERM Convert spectral net flux to the thermal grid.
   %
   % Qsi - timestep total incoming solar radiation flux
   % I0 - total incoming solar radiation flux used in two-stream solve
   % albedo - timestep (actual) albedo
   % Qnet - interface net solar radiation flux on the spectral mesh
   % dz_spect - subsurface spectral mesh layer thicknesses
   % dz_therm - subsurface thermal mesh layer thicknesses

   % Convert the net-flux profile to absorbed flux over each spectral CV.
   dQnet_spect = Qnet(1:end - 1) - Qnet(2:end);

   % Aggregate the absorbed spectral flux onto the thermal control volumes.
   % This current collapse assumes the thermal/spectral spacing ratio is fixed
   % and uniform. If the spectral grid becomes nonuniform, revisit this mapping
   % instead of only changing icemodel.radiation.bulk_extinction_coefficients.
   n_spect_per_therm = dz_therm(1) / dz_spect(1);
   dQnet_therm = transpose(sum(reshape(dQnet_spect, n_spect_per_therm, []), 1));
   dQnet_therm = [dQnet_therm; zeros((sum(dz_therm) - sum(dz_spect)) ...
      / dz_therm(1), 1)];

   % Preserve the established chi rule used by the SEB coupling.
   if albedo > 0.65
      chi = 0.9;
   else
      chi = dQnet_therm(1) / sum(dQnet_therm);
   end

   % Convert absorbed flux to a volumetric thermal source term. The Q profile
   % is scaled by the prototype spectrum I0, so rescale it back to the actual
   % timestep shortwave forcing Qsi here.
   Sc = (1.0 - chi) * Qsi / I0 * -dQnet_therm ./ dz_therm;
end

%--------------------------------------------------------------------------
% Notes
%
% The two-stream solution gives the total upflux (X) and downflux (Y) from
% the respective boundary to each layer interface. Qnet (XYnet in Schlatter) is
% the net flux formed from those interface fluxes, and dQ = d/dz(Qnet) is the
% absorbed flux.
%
%   Y↓(1)  X↑(1)    Xnet = X↑(2)-X↑(1), Ynet = Y↓(1)-Y↓(2)
%   ____________
%   ____________    Qnet = Xnet - Ynet = (X↑(2)-X↑(1))-(Y↓(1)-Y↓(2))
%   Y↓(2)  X↑(2)
%
% Rescaling to forcing Qsi:
%
% For the prototype spectrum, the depth-integrated absorbed flux is:
% (1) -sum(dQ) = (1 - albedo) * I0
%
% but the absorbed source term for the model timestep should satisfy:
% (2) -sum(dQ) = (1 - albedo) * Qsi
%
% so multiply (1) by Qsi / I0:
% (3) (1 - albedo) * I0 * Qsi / I0 = (1 - albedo) * Qsi
%
% and rescale the absorbed spectral profile by Qsi / I0 when constructing Sc:
%
% (4) Sc = Qsi / I0 .* -dQ ./ dz
%
% if I0 were the actual incoming solar at each timestep, it would be:
% (5) Sc = -dQ ./ dz
%
% To assign a portion of the incoming radiation, chi, to surface heating, and
% the rest to subsurface layers, Sc need to be adjusted for the surface portion,
% hence the implementation:
%
% Sc = (1.0 - chi) * Qsi / I0 * -dQ_therm ./ dz_therm;
%
% The surface portion is allocated in the SEB. Note that Qnet and dQnet are not
% adjusted for chi in this routine.
%--------------------------------------------------------------------------
