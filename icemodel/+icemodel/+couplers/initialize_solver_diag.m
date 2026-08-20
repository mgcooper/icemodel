function diag = initialize_solver_diag()
   %INITIALIZE_SOLVER_DIAG Return the sentinel solver-diagnostics struct.
   %
   %  diag = icemodel.couplers.initialize_solver_diag()
   %
   % One fixed-schema scalar struct carries every solver observability
   % value from both surface-column couplers to the driver. The `ok`
   % control flags stay plain coupler returns; this struct exists so the
   % diagnostics can change without touching coupler signatures.
   %
   % The driver holds one sentinel copy per forcing step and overwrites
   % it only at the accepted-substep checkpoint, so a forcing step whose
   % only content was forced advances reports these sentinel values,
   % never a rejected solve's values. The convergence flags and the inner
   % iteration count are output-visible in standard-profile channels;
   % cpl_iters, cpl_res, and seb_res are emitted only through the
   % diagnostic-profile suffix. The sentinels are therefore fixed here:
   % converged flags false, iteration counts and residuals NaN or zero.
   %
   % Fields
   %  ok_seb, ok_ieb, ok_cpl - copies of the control flags [logical]
   %  n_iters      - inner Picard iterations, final sweep
   %  cpl_iters    - outer iterations used, all phases together
   %  cpl_phase    - 0 sentinel, 1 primary accepted, 2 conservative
   %  cpl_recovered - the conservative phase produced the returned state
   %  cpl_res      - final outer residual T_sfc - Ts_old [K]
   %  seb_res      - final SEB residual magnitude [W m-2]
   %  cpl_res_hist - signed outer-residual ring, newest last (16 x 1);
   %                 rides the debug dumps only, never an output channel
   %
   % See also: icemodel.couplers.solve_surface_column_robin,
   %  icemodel.couplers.solve_surface_column_dirichlet
   %
   %#codegen

   diag = struct( ...
      'ok_seb', false, ...
      'ok_ieb', false, ...
      'ok_cpl', false, ...
      'n_iters', nan, ...
      'cpl_iters', 0.0, ...
      'cpl_phase', 0.0, ...
      'cpl_recovered', false, ...
      'cpl_res', nan, ...
      'seb_res', nan, ...
      'cpl_res_hist', nan(16, 1));
end
