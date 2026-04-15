function [De, br_coefs, liqflag, ro_air_Lv, ea_atm, ro_sfc] = ...
      initialize_surface_state(f_liq, f_ice, wspd, tair, rh, opts)
   %INITIALIZE_SURFACE_STATE Initialize the surface state for the SEB solver.
   %
   %  [De, br_coefs, liqflag, ro_air_Lv, ea_atm] = ...
   %     icemodel.surface.initialize_surface_state(f_liq, f_ice, wspd, tair, rh, opts)
   %  [De, br_coefs, liqflag, ro_air_Lv, ea_atm, ro_sfc] = ...
   %     icemodel.surface.initialize_surface_state(f_liq, f_ice, wspd, tair, rh, opts)
   %
   %  Initializes all surface-state quantities derived from forcing data
   %  rather than the ice-column profile:
   %
   %  De         — pre-computed bulk-Richardson aerodynamic exchange coefficient
   %               for every forcing step [m s-1]
   %  br_coefs   — pre-computed bulk-Richardson stability-factor coefficients;
   %               scalar constants [gamma, S2, S3]
   %  liqflag    — bootstrap surface phase flag; true if the initial surface
   %               layer contains liquid water
   %  ro_air_Lv  — initial ro_air × latent heat [J m-3], matching liqflag
   %               (roLv if liqflag, roLs otherwise); updated each substep
   %               by icemodel.timestepping.updatesubstep
   %  ea_atm     — pre-computed atmospheric vapor pressure for every forcing
   %               step [Pa], computed using tair > Tf as the phase selector
   %               (WMO convention: RH data are referenced to liquid water
   %               for T > 0 °C and to ice for T < 0 °C)
   %  ro_sfc     — initial surface bulk density [kg m-3]; updated each substep
   %               by icemodel.timestepping.updatesubstep
   %
   %  Centralizing these initializations ensures that:
   %   (a) the exchange-coefficient loop over the full wspd vector runs once;
   %   (b) ea_atm uses the physically correct tair-based phase flag;
   %   (c) icemodel, skinmodel, and the test fixtures share one canonical
   %       initialization path.
   %
   %  Inputs:
   %    f_liq   — initial volumetric liquid-water fraction profile [JJ × 1]
   %    f_ice   — initial volumetric ice fraction profile [JJ × 1]
   %    wspd    — wind-speed forcing for all timesteps [N × 1], m s-1
   %    tair    — air-temperature forcing for all timesteps [N × 1], K
   %    rh      — relative humidity forcing for all timesteps [N × 1], %
   %    opts    — model options struct (uses z0_bulk, z_tair, z_wind)
   %
   % See also:
   %   icemodel.surface.turbulence.bulk_richardson.exchange_coefficients,
   %   icemodel.surface.atmospheric_vapor_pressure,
   %   icemodel.timestepping.updatesubstep
   %
   %#codegen

   persistent Tf roLs roLv f_liq_phase_switch_threshold
   if isempty(Tf)
      [Tf, roLs, roLv] = icemodel.physicalConstant('Tf', 'roLs', 'roLv');
      f_liq_phase_switch_threshold = icemodel.parameterLookup( ...
         'f_liq_phase_switch_threshold');
   end

   % Pre-compute bulk-Richardson exchange coefficients for every forcing step.
   [De, br_coefs] = ...
      icemodel.surface.turbulence.bulk_richardson.exchange_coefficients( ...
      wspd, opts.z0_bulk, opts.z_tair, opts.z_wind);

   % Bootstrap surface phase flag from the initial surface liquid fraction.
   liqflag = f_liq(1) > f_liq_phase_switch_threshold;
   if liqflag
      ro_air_Lv = roLv;
   else
      ro_air_Lv = roLs;
   end

   % Pre-compute atmospheric vapor pressure for every forcing step.
   ea_atm = icemodel.surface.atmospheric_vapor_pressure(tair, rh);

   % Bootstrap surface bulk density from the initial surface layer.
   ro_sfc = icemodel.surface.surface_bulk_density(f_ice(1), f_liq(1));
end
