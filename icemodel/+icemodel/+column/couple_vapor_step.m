function [f_ice, f_liq, d_sbl_err, ledger] = couple_vapor_step( ...
      T, f_ice, f_liq, d_sbl_err, dz, delz, fn, dt, f_ice_min, ...
      f_res_por, ledger, use_mass_budget)
   %COUPLE_VAPOR_STEP Move vapor between the cells and apply it, once.
   %
   %  [f_ice, f_liq, d_sbl_err, ledger] = ...
   %     icemodel.column.couple_vapor_step( ...
   %     T, f_ice, f_liq, d_sbl_err, dz, delz, fn, dt, f_ice_min, ...
   %     f_res_por, ledger, use_mass_budget)
   %
   % Call this once per accepted substep, before the surface mass balance.
   % Interior transport is a separate step from the surface exchange because
   % the two conserve different things: Fick's law fixes the mass here, while
   % the surface energy balance fixes the energy there.
   %
   % The whole coupled interior path lives behind this one call, so the driver
   % holds no vapor intermediates. It evaluates the accepted-state node
   % quantities. It moves vapor across the interior faces, applies the
   % arriving mass under the three per-cell limits, and records what the
   % ledger needs.
   %
   % D_SBL_ERR threads in and out. The applier records what a limited cell
   % could not take. This adds that to the running per-cell total, so the
   % increment accumulates the way every other d_* increment does.
   %
   % This evaluates the node quantities at the accepted substep state. The
   % grain-growth call in icemodel.m evaluates its own at the end-of-step
   % state, where T differs. The two are not duplicates and must not be
   % merged.
   %
   % Inputs
   %   T              - Column temperature at the accepted substep [K].
   %   f_ice, f_liq   - Column phase fractions [-].
   %   d_sbl_err      - Running per-cell unapplied vapor record [-], the
   %                    ice-fraction equivalent on the energy basis.
   %   dz, delz, fn   - Control-volume thickness, node spacing, and interface
   %                    weights.
   %   dt             - Substep length [s].
   %   f_ice_min      - Minimum retained ice fraction [-].
   %   f_res_por      - Residual liquid-water fraction per pore volume [-].
   %   ledger         - Forcing-step mass and energy ledger.
   %   use_mass_budget - True when the diagnostic ledger is being built.
   %
   % Outputs
   %   f_ice, f_liq   - Phase fractions after the interior exchange.
   %   d_sbl_err      - The running record with this step's shortfall added.
   %   ledger         - The ledger with the redistribution energy accumulated.
   %
   % See also: icemodel.column.couple_vapor_transport,
   %  icemodel.column.apply_vapor_transport,
   %  icemodel.column.accumulate_redistribution_budget
   %
   %#codegen

   % Node quantities at the accepted state, evaluated once and reused by the
   % transport below. That is one exponential and one power per accepted
   % substep, not one per solver iteration.
   [ro_vap, De] = icemodel.column.accepted_vapor_quantities(T, f_liq);

   % Fick's law across the interior faces. Both boundaries are closed, so
   % this redistributes and creates nothing.
   d_vap = icemodel.column.couple_vapor_transport( ...
      ro_vap, De, dz, delz, fn, dt);

   % Storage before the exchange, for the redistribution energy below.
   if use_mass_budget
      [solid_r, liquid_r] = icemodel.column.integrate_column_budget( ...
         T, f_ice, f_liq, dz);
   end

   [f_ice, f_liq, d_sbl_err_cpl] = icemodel.column.apply_vapor_transport( ...
      f_ice, f_liq, d_vap, f_ice_min, f_res_por);

   % Carry the unapplied amount into the running record, so nothing the
   % interior could not take goes unrecorded.
   d_sbl_err = d_sbl_err + d_sbl_err_cpl;

   % Interior transport conserves mass, but it moves that mass between cells
   % of different phase. The latent heats differ, so the energy-weighted
   % storage moves with no potential input. Record that so the vapor closure
   % identity still balances.
   if use_mass_budget
      ledger = icemodel.column.accumulate_redistribution_budget( ...
         ledger, solid_r, liquid_r, T, f_ice, f_liq, dz);
   end
end
