function [T, f_ice, f_liq, d_liq, d_evp, d_drn, x_err, lcflag] = ICEMF( ...
      T, f_ice, f_liq, ro_ice, ro_liq, cv_ice, cv_liq, Lf, Ls, Lv, Tf, ...
      TL, fcp, xf_liq, Sc, Sp, JJ, f_min, dz_therm, dt_new, Qe, ro_iwe, ...
      d_liq, d_drn, d_evp, flmin, liqresid)
   %ICEMF compute ice melt-freeze and combine layers if necessary
   %
   %#codegen

   % Compute d_liq first. xf_liq is the value of f_liq returned by this function
   % at the end of the last step. The only process that affects f_liq between
   % its assignment to xf_liq and here is melt/freeze in ICEENBAL. Therefore,
   % d_liq represents the change in f_liq due to melt/freeze.
   % Evaporation/condensation is accounted for below. d_drn gets updated if
   % layers are combined
   % liqsum(1) = sum(f_liq.*dz_therm);
   % liqsum(2) = nan;

   d_liq = d_liq + f_liq - xf_liq; % d_liq < 0 means refreezing, d_liq > 0 melt
   xf_liq = f_liq;
   % xf_ice = f_ice;

   % want:
   % d_liq > 0 = melt (should be correct below)
   % d_evp > 0 = cond (should be correct below)
   % d_drn > 0 = runoff due to extra condensation (should be correct below)

   % Update the lower melt-zone boundary
   if T(1) > TL
      
      % Ensure f_liq is never reduced below residual water for melting nodes
      f_liq_min = (f_liq(1) + f_ice(1) * ro_iwe) * flmin;
      liqresid = max(liqresid, f_liq_min / (1 - f_ice(1)));

      % If f_res is larger than f_liq_min, ICESUBL will never reduce f_liq(1)
      % below f_liq_min. Thus, below checks if f_res is ever smaller than
      % f_liq_min. If it is, then there's a problem 
      
      % Check if f_res is ever smaller than f_liq_min. If it can be proven it
      % will never trigger, delete this entire if block and use liq_resid
      % without any max() statement.
      % 
      % if liq_resid * (1 - f_ice(1)) < f_liq_min
      %    liq_resid = max(0.07, f_liq_min / (1 - f_ice(1)));
      % end
   end

   % evaporation/condensation/sublimation
   [f_ice, f_liq, d_drn, x_err] = ICESUBL(f_ice, f_liq, d_drn, ...
      ro_ice, ro_liq, Ls, Lv, dz_therm, dt_new, Qe, f_min, liqresid);

   % % budget evap / subl
   d_evp = d_evp + f_liq - xf_liq;
   % d_sub = f_ice - xf_ice;
   % if evap, d_evp < 0
   % if cond, d_evp > 0
   % if subl, d_sub < 0

   % Update the top layer temperature, if the node is in the melt zone
   if T(1) > TL % f_liq(1) > f_liq_min
      T(1) = Tf - sqrt( (f_liq(1) + f_ice(1) * ro_iwe) / f_liq(1) - 1.0) / fcp;
   end

   % Below here:
   % d_XXX = f_liq - xf_liq - d_evp
   % where f_liq on the rhs is f_liq_new
   
   % combine layers if any layer is <f_min, or if this step's sublimation
   % would reduce any layer to <f_min (predict the need to combine next step)
   lyrmrg = f_ice <= f_min | (f_ice + Qe/(Ls*ro_ice)*dt_new/dz_therm) <= f_min;

   % lyrmrg is updated in the loop, keep lcflag to know which layer was combined
   lcflag = lyrmrg;

   % if lcflag is all false, return
   if any(lyrmrg)
      % This is a counter that keeps the ji index on track with the original
      % 1:JJ column, since layers are removed within the loop
      ii = 0;

      for j = 1:JJ
         ji = j + ii;

         if lyrmrg(ji)
            % Combine layers
            [j1,j2] = LAYERINDS(ji,f_ice);

            [f_liq(j2), f_ice(j2), T(j2), Sc(j2), Sp(j2), d_drn, d_liq] = ...
               COMBINEHEAT(f_liq, f_ice, T, Sc, Sp, Tf, TL, cv_ice, cv_liq, ...
               Lf, ro_ice, ro_liq, fcp, j1, j2, d_drn, d_liq, dz_therm);

            % Remove the combined layers and add new layers to the bottom.
            T = vertcat(T, T(end)); T(j1) = [];
            Sc = vertcat(Sc, Sc(end)); Sc(j1) = [];
            Sp = vertcat(Sp, Sp(end)); Sp(j1) = [];
            f_ice = vertcat(f_ice, f_ice(end)); f_ice(j1) = [];
            f_liq = vertcat(f_liq, f_liq(end)); f_liq(j1) = [];
            lyrmrg = vertcat(lyrmrg, lyrmrg(end)); lyrmrg(j1) = [];
            
            ii = ii - 1;
         end
      end

      % try this - it's the original one adjusted for evap
      d_drn = xf_liq - f_liq + d_evp;

      % liqsum(2) = sum(f_liq.*dz_therm);
      % % compute d_drn due to layer combination:
      % d_com = f_liq - xf_liq - d_evp;
      %
      % % allocate drained water to d_drn
      % d_drn(d_com < 0) = d_drn(d_com < 0) - d_com(d_com < 0);

      % otherwise, any change in f_liq is accounted for on the next timestep

      % think this is wrong
      % d_drn = d_drn + f_liq - xf_liq - d_evp; % = d_drn + d_com

      % original:
      % d_drn = d_drn + xf_liq - f_liq;

      % test 8
      % d_drn = d_drn + f_liq - xf_liq;

      % d_drn = d_drn + f_liq - xf_liq - d_evp;
   end
end

% % f_liq budget:
% 1 d_liq = f_liq - xf_liq (melt / freeze - if d_liq < 0, it means freeze)
% 2 xf_liq = f_liq
% 3 f_liq = f_liq + d_evap (evap / cond in ICESUBL)
% 4 d_drn = d_drn + x_evap (extra condensation that runs off in ICESUBL)
% 5 f_liq = f_liq + d_comb (extra water that runs off in COMBINEHEAT)
% 6 d_drn = d_drn + xf_liq - f_liq
%
% Combining terms, step 6 can be written:
% d_drn = x_evap + f_liq - (f_liq _ d_evap + d_comb)
% d_drn = x_evap - d_evap - d_comb
% d_drn = x_evap - d_evap - d_comb
%
% if x_evap is non-zero, it is guaranteed to be positive, and it means
% condensation occurred, therefore d_evap is positive, and d_evap


% % f_liq budget:
% 1 d_liq = f_liq - xf_liq
% 2 xf_liq = f_liq
% 3 f_liq = f_liq + d_evap (evap / cond in ICESUBL)
% 4 d_drn = d_drn + x_evap (extra condensation that runs off in ICESUBL)
% 5 f_liq = f_liq + d_comb (extra water that runs off in COMBINEHEAT)
% 6 d_drn = d_drn + xf_liq - f_liq
%
% Combining terms, step 6 can be written:
% d_drn = x_evap + f_liq - (f_liq _ d_evap + d_comb)
% d_drn = x_evap - d_evap - d_comb
% d_drn = x_evap - d_evap - d_comb
%
% if x_evap is non-zero, it is guaranteed to be positive, and it means
% condensation occurred, therefore d_evap is positive, and d_evap

% assert(f_liq(1) + f_ice(1) <= 1.0)
% printf(f_liq(1), 10)
% printf(f_ice(1), 10)
% printf(T(1), 10)
% printf(f_liq_min(1), 10)

% xf_liq is from the prior met-step, so each sub-step d_liq is overridden,
% and the total d_liq is relative to each met-step, as it should be,
% whereas the other values in here are updated at each substep. However,
% what happens when a layr is removed within a sub-stepping?
