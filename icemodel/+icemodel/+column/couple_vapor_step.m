function [f_ice, f_liq, d_vap_faces, vapor_solid, vapor_liquid, ledger] = ...
      couple_vapor_step(T, f_ice, f_liq, f_ice_solve, f_liq_solve, ...
      d_vap_faces, vapor_solid, vapor_liquid, dz, delz, fn, dt, ...
      f_ice_min, f_res_por, ledger, use_mass_budget)
   %COUPLE_VAPOR_STEP Move vapor between the cells and apply it, once.
   %
   %  [f_ice, f_liq, d_vap_faces, vapor_solid, vapor_liquid, ledger] = ...
   %     icemodel.column.couple_vapor_step( ...
   %     T, f_ice, f_liq, f_ice_solve, f_liq_solve, d_vap_faces, ...
   %     vapor_solid, vapor_liquid, dz, delz, fn, dt, f_ice_min, ...
   %     f_res_por, ledger, use_mass_budget)
   %
   % Call this once per accepted substep, after the surface mass balance and
   % its vapor budget close. The position matters three ways. The surface
   % path's d_liq increment spans the checkpoint to the surface budget call,
   % so transport before it would be read as melt or refreezing. The vapor
   % budget's storage baseline spans the surface exchange alone, so
   % transport before it would be scored as surface exchange. And the
   % shortfall a limited cell records here goes to the redistribution
   % accounting below, never into the surface closure identity.
   %
   % Interior transport is a separate step from the surface exchange because
   % the two conserve different things: Fick's law fixes the mass here, while
   % the surface energy balance fixes the energy there.
   %
   % The whole coupled interior path lives behind this one call. It moves
   % vapor across the interior faces, applies the arriving mass under the
   % three per-cell limits, and records what the ledger needs.
   %
   % F_ICE_SOLVE and F_LIQ_SOLVE are the fractions the solve converged on,
   % before the surface exchange. The face quantities AND the per-cell
   % phase decisions evaluate at that state, so the mass this moves, the
   % phase it moves between, and the latent heats it implies are the
   % conjugates of the energy the solve transported. The exchange can flip
   % the top cell's wet/dry class; deciding from the post-exchange state
   % would let the energy carry one phase's latent heat while the applier
   % spends the other's. The per-cell limits still clamp against the
   % current F_ICE and F_LIQ, because a cell can only give what it holds
   % now.
   %
   % D_VAP_FACES threads in and out: the gross water-equivalent depth each
   % face exchanged this forcing step. This adds the substep's interior face
   % magnitudes, so grain growth can consume the fluxes the column actually
   % transported rather than a snapshot recomputed at end-of-step state.
   % Face 1 belongs to the driver, which adds the realized surface exchange;
   % the transport's own face 1 is closed and adds zero.
   %
   % VAPOR_SOLID and VAPOR_LIQUID thread in and out: the substep's
   % storage-change context [m w.e.] that the follow-on remesh budget uses
   % for the endpoint-storage gross channels. This adds the transport's
   % per-phase increments, so a transport-only substep cannot report a zero
   % endpoint gross while its per-phase endpoints move.
   %
   % This evaluates the node quantities at the solve state once per
   % accepted substep, one exponential and one power, not one per solver
   % iteration; the grain-growth call needs none at all because it
   % consumes the accumulation this builds.
   %
   % Inputs
   %   T              - Column temperature at the accepted substep [K].
   %   f_ice, f_liq   - Column phase fractions, after the surface
   %                    exchange [-].
   %   f_ice_solve    - Ice fraction the solve converged on [-].
   %   f_liq_solve    - Liquid fraction the solve converged on [-].
   %   d_vap_faces    - Running gross face exchange [m w.e.] (JJ+1 x 1).
   %   vapor_solid    - Substep solid storage-change context [m w.e.].
   %   vapor_liquid   - Substep liquid storage-change context [m w.e.].
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
   %   d_vap_faces    - The running record with this substep's interior face
   %                    magnitudes added.
   %   vapor_solid,   - The storage-change context with the transport's
   %   vapor_liquid     per-phase increments added.
   %   ledger         - The ledger with the per-phase redistribution
   %                    increments and the interior shortfall accumulated.
   %
   % See also: icemodel.column.couple_vapor_transport,
   %  icemodel.column.apply_vapor_transport,
   %  icemodel.column.accumulate_redistribution_budget
   %
   %#codegen

   persistent ro_liq
   if isempty(ro_liq)
      ro_liq = icemodel.physicalConstant('ro_liq');
   end

   % Node quantities at the state the solve converged on, evaluated once
   % and reused by the transport below. Not the post-exchange fractions:
   % the solve's face energy came from this state, and conjugacy needs the
   % mass to come from the same one.
   [ro_vap, De] = icemodel.column.accepted_vapor_quantities(T, f_liq_solve);

   % One phase decision per cell, at the solve state, shared by the
   % applier below so the phase the mass lands in matches the latent heat
   % the solve's face energy carried.
   [wet, f_res] = icemodel.column.vapor_exchange_is_wet( ...
      f_ice_solve, f_liq_solve, f_res_por);

   % Fick's law across the interior faces. Both boundaries are closed, so
   % this redistributes and creates nothing.
   [d_vap, ~, U_vap_faces] = icemodel.column.couple_vapor_transport( ...
      ro_vap, De, dz, delz, fn, dt);

   % Add this substep's face magnitudes to the step accumulation, as
   % water-equivalent depth. Gross, because grain growth scales with the
   % flux magnitude and reversing substeps must add rather than cancel.
   d_vap_faces = d_vap_faces + abs(U_vap_faces) * dt / ro_liq;

   % Storage before the exchange, for the redistribution increments below.
   if use_mass_budget
      [solid_r, liquid_r] = icemodel.column.integrate_column_budget( ...
         T, f_ice, f_liq, dz);
   end

   [f_ice, f_liq, d_sbl_err_cpl] = icemodel.column.apply_vapor_transport( ...
      f_ice, f_liq, d_vap, f_ice_min, wet, f_res);

   % Interior transport conserves mass, but it moves that mass between
   % cells of different phase, so the solid and liquid storage totals move
   % in opposite directions while their sum holds. Record the per-phase
   % increments and the shortfall the per-cell limits rejected in the
   % redistribution channels, and add the increments to the substep
   % storage-change context the remesh budget reads. The surface closure
   % identity stays surface-only; the storage closures consume the
   % per-phase increments.
   if use_mass_budget
      [solid_after, liquid_after] = icemodel.column.integrate_column_budget( ...
         T, f_ice, f_liq, dz);
      vapor_solid = vapor_solid + (solid_after - solid_r);
      vapor_liquid = vapor_liquid + (liquid_after - liquid_r);
      ledger = icemodel.column.accumulate_redistribution_budget( ...
         ledger, solid_r, liquid_r, T, f_ice, f_liq, dz, d_sbl_err_cpl);
   end
end
