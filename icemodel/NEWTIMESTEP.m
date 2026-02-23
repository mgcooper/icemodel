function [dt_sum, subfail, ok_seb, ok_ieb, d_liq, d_evp, d_lyr] = ...
      NEWTIMESTEP(f_liq, opts)
   %NEWTIMESTEP initialize new timestep
   %
   %#codegen

   ok_seb = false;         % assume seb failed
   ok_ieb = false;         % assume ice-eb sub-step failed
   d_liq = 0.0 * f_liq;    % reset the change in liq water content
   d_evp = 0.0 * f_liq;    % reset the evaporation change in water content
   d_lyr = 0.0 * f_liq;    % reset the layer change
   dt_sum = 0.0;
   subfail = 0;            % keep track of failed substeps

   % For non-Dirichlet bc, set ok_seb true so dt control advances (see NEXTSTEP)
   if opts.bc_type > 1
      ok_seb = true;
   end
end
