function diag = update_solver_diag(cpl_recovery_active, diag)
   %UPDATE_SOLVER_DIAG Record the accepted substep solver result.
   %
   %  diag = icemodel.couplers.update_solver_diag( ...
   %     cpl_recovery_active, diag)
   %
   % CPL_RECOVERY_ACTIVE is true when the accepted solve used recovery mode.
   % DIAG contains the forcing-step fields and the latest attempt in
   % DIAG.SUBSTEP. Copy the attempt only when all three solver checks pass.
   % Keep the forcing-step counters already in DIAG.
   %
   % See also: icemodel.couplers.initialize_solver_diag,
   %  icemodel.timestepping.acceptsubstep
   %
   if diag.substep.ok_seb && diag.substep.ok_ieb && diag.substep.ok_cpl
      diag = icemodel.helpers.copyFields(diag, diag.substep, true);
      diag.cpl_recovery_count = ...
         diag.cpl_recovery_count + double(cpl_recovery_active);
   end
end
