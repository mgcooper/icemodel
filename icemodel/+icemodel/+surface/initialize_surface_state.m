function [ea_atm, ro_atm, cv_atm, nu_air, H_h, De_e, br_coefs] ...
      = initialize_surface_state(opts, tair, wspd, rh, psfc)
   %INITIALIZE_SURFACE_STATE Precompute forcing-derived surface arrays.
   %
   %  [ea_atm, ro_atm, cv_atm, nu_air, H_h, De_e, br_coefs] = ...
   %     icemodel.surface.initialize_surface_state(opts, tair, wspd, rh, psfc)
   %
   %  Precomputes all surface-state quantities that are pure functions
   %  of the meteorological forcing and model geometry. These arrays
   %  are indexed by metstep in the main flow and remain constant for
   %  the duration of the simulation.
   %
   %  ea_atm  — atmospheric vapor pressure [Pa], one per forcing step
   %  ro_atm  — moist-air density [kg m-3], one per forcing step
   %  cv_atm  — volumetric heat capacity of moist air [J m-3 K-1]
   %  nu_air  — kinematic viscosity of air [m2 s-1]
   %  H_h     — sensible heat transport prefactor [W m-2 K-1]
   %            = cv_atm .* De_h
   %  De_e    — latent exchange coefficient [m s-1 Pa-1]
   %            = De_h .* epsilon ./ psfc
   %  br_coefs — bulk-Richardson stability coefficients [gamma S2 S3]
   %
   %  Surface running state (liqflag, ro_sfc, hv_atm, H_e, f_res_por)
   %  is NOT computed here. Those quantities depend on the evolving
   %  column state and current forcing step, so they are derived at
   %  each substep entry by icemodel.surface.update_surface_state.
   %
   %  Inputs:
   %    opts  — model options struct (uses z0_bulk, z_tair, z_wind)
   %    tair  — air-temperature forcing [N x 1], K
   %    wspd  — wind-speed forcing [N x 1], m s-1
   %    rh    — relative humidity forcing [N x 1], %
   %    psfc  — surface pressure forcing [N x 1], Pa
   %
   % See also:
   %   icemodel.surface.update_surface_state,
   %   icemodel.surface.turbulence.bulk_richardson.exchange_coefficients,
   %   icemodel.surface.atmospheric_vapor_pressure
   %
   %#codegen

   persistent cp_air epsilon
   if isempty(cp_air)
      [cp_air, epsilon] = icemodel.physicalConstant( ...
         'cp_air', 'epsilon');
   end

   % Pre-compute bulk-Richardson exchange coefficients for every
   % forcing step. De_h is the aerodynamic exchange coefficient for
   % heat [m s-1], used here only to derive H_h and De_e. br_coefs
   % are the Louis, 1979 stability coefficients.
   [De_h, br_coefs] = ...
      icemodel.surface.turbulence.bulk_richardson. ...
      exchange_coefficients( ...
      wspd, opts.z0_bulk, opts.z_tair, opts.z_wind);

   % Pre-compute atmospheric vapor pressure for every forcing step.
   ea_atm = icemodel.surface.atmospheric_vapor_pressure(tair, rh);

   % Moist-air density and kinematic viscosity for every forcing step.
   ro_atm = icemodel.vapor.moist_air_density(psfc, ea_atm, tair);
   nu_air = icemodel.kernels.air_kinematic_viscosity(tair, ro_atm);

   % Volumetric heat capacity of moist air [J m-3 K-1]
   cv_atm = ro_atm * cp_air;

   % Sensible heat transport prefactor [W m-2 K-1]
   H_h = cv_atm .* De_h;

   % Latent exchange coefficient [m s-1 Pa-1]
   De_e = De_h .* epsilon ./ psfc;
end
