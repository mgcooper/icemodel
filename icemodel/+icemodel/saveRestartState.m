function saveRestartState(opts, simyear, T_ice, f_ice, f_liq, T_sfc, r_eff)
   %SAVERESTARTSTATE Save the year-boundary state needed for a restart.
   %
   %  icemodel.saveRestartState(opts, simyear, T_ice, f_ice, f_liq, T_sfc, ...
   %     r_eff)
   %
   %  The restart file preserves variables that evolve through the governing
   %  equations. Rebuilding them requires a rerun of the prior simulation:
   %
   %    T_ice  - ice column temperature [K]
   %    f_ice  - volumetric ice fraction [-]
   %    f_liq  - volumetric liquid fraction [-]
   %    T_sfc  - surface temperature [K] (SEB boundary condition / prior)
   %    r_eff  - effective grain radius [m] (vapor mass transfer)
   %
   %  Each field takes the name of its model variable.
   %
   %  checksubstep selects retry settings for each substep. acceptsubstep
   %  restores the primary settings. The restart file does not store solver
   %  settings.
   %
   %  Surface running state (liqflag, ro_sfc, hv_atm, H_e, f_res_por) is
   %  NOT saved. Those quantities are derived at substep entry by
   %  icemodel.surface.update_surface_state from the current column state
   %  and forcing step, so they are always consistent and need not be
   %  preserved across restart boundaries.
   %
   %  The restart file also stores metadata (casename, smbmodel, opts,
   %  timestamp) for provenance and validation by loadRestartState.
   %
   %  See also: icemodel.loadRestartState,
   %     icemodel.column.initialize_column_state,
   %     icemodel.surface.update_surface_state

   restart = struct();
   restart.simyear = simyear;

   % Prognostic state
   restart.T_ice = T_ice;
   restart.f_ice = f_ice;
   restart.f_liq = f_liq;
   restart.T_sfc = T_sfc;
   restart.r_eff = r_eff;

   % Metadata
   restart.casename = string(opts.casename);
   restart.smbmodel = string(opts.smbmodel);
   restart.opts = opts;
   restart.timestamp_utc = datetime('now', 'TimeZone', 'UTC');

   filepath = icemodel.restartfile(opts, simyear);
   backupfile(filepath, opts.backupflag);
   save(filepath, 'restart');
end
