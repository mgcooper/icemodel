function [T, f_ice, f_liq, d_liq, d_evp, d_lyr, x_err, lcflag] = ICEMF( ...
      T, f_ice, f_liq, ro_ice, ro_liq, cv_ice, cv_liq, Lf, Ls, Lv, Tf, ...
      TL, fcp, xf_liq, Sc, Sp, JJ, f_min, dz_therm, d_pevp, d_liq, ...
      d_evp, d_lyr, f_ell_min, f_liq_res)
   %ICEMF Budget ice melt/freeze, evap/cond, and remesh melted layers.
   %
   % Arguments
   %  T - ice temperature
   %  f_ice - ice fraction
   %  f_liq - liquid fraction
   %  d_liq - the total change in liquid water fraction over the full timestep.
   %  xf_liq - the prior value of f_liq returned by this function.
   %  d_evp - the change in liquid water fraction due to evaporation.
   %  d_con - the change in liquid water fraction due to condensation.
   %  d_pevp - potential vapor-driven change in top-layer liquid fraction.
   %  x_err - latent heat flux which exceeds the top layer liq/ice content.
   %  x_err is only used for debugging, it is not included in any mass budgets.
   %
   % See also:
   %
   %#codegen

   % Compute the surface mass balance terms
   [T, d_liq, d_evp, ~, f_ice, f_liq, x_err] = SMB(T, d_liq, d_evp, 0, ...
      f_ice, f_liq, xf_liq, TL, ro_ice, ro_liq, f_ell_min, f_liq_res, ...
      Ls, Lv, d_pevp, f_min, fcp, Tf);

   % Combine layers if ice fraction is <f_min, or if this step's sublimation
   % would reduce it below <f_min (i.e., predict the need to combine next step)
   lyrmrg = f_ice <= f_min | ...
      (f_ice + d_pevp * (Lv * ro_liq) / (Ls * ro_ice)) <= f_min;

   % lyrmrg is updated in the loop, lcflag indicates which layer(s) combined
   lcflag = lyrmrg;

   if any(lyrmrg)
      % This is a counter that keeps the ji index on track with the original
      % 1:JJ column, since layers are removed within the loop
      ii = 0;

      for j = 1:JJ
         ji = j + ii;

         if lyrmrg(ji)
            % Combine layers
            [j1,j2] = LAYERINDS(ji,f_ice);

            [T(j2), f_ice(j2), f_liq(j2), Sc(j2), Sp(j2), d_lyr] = ...
               COMBINEHEAT(T, f_ice, f_liq, Sc, Sp, Tf, TL, cv_ice, cv_liq, ...
               Lf, ro_ice, ro_liq, fcp, j1, j2, d_lyr, dz_therm);

            % Remove the combined layers and add new layers to the bottom.
            T = vertcat(T, T(end)); T(j1) = []; %#ok<*AGROW>
            Sc = vertcat(Sc, Sc(end)); Sc(j1) = [];
            Sp = vertcat(Sp, Sp(end)); Sp(j1) = [];
            f_ice = vertcat(f_ice, f_ice(end)); f_ice(j1) = [];
            f_liq = vertcat(f_liq, f_liq(end)); f_liq(j1) = [];
            lyrmrg = vertcat(lyrmrg, lyrmrg(end)); lyrmrg(j1) = [];

            ii = ii - 1;
         end
      end
   end
end

function [T, d_liq, d_evp, d_con, f_ice, f_liq, x_err] = SMB(T, d_liq, d_evp, ...
      d_con, f_ice, f_liq, xf_liq, TL, ro_ice, ro_liq, f_ell_min, f_liq_res, ...
      Ls, Lv, d_pevp, f_min, fcp, Tf)
   %SMB (sub)surface mass budgets
   %
   % Positive values of d_liq indicate increasing water fraction:
   % d_liq > 0 = melt.
   % d_evp > 0 = condensation.
   % d_con > 0 = condensation which exceeds top layer porosity.
   % d_liq < 0 = freeze.
   % d_evp < 0 = evaporation.

   % Compute delta f_liq. Note, the only process which affects f_liq between
   % d_liq assignments here is phase change (ICEENBAL). d_liq < 0 means
   % refreezing, d_liq > 0 melt.
   d_liq = d_liq + f_liq - xf_liq;

   % Reset past values for budgeting evap/condensation.
   xf_liq = f_liq;
   xf_ice = f_ice;

   % Update the residual unfrozen water fraction defined by the phase fraction
   % characteristic function at the lower temperature boundary. This ensures
   % f_liq is not reduced below residual water for melting nodes.
   if T(1) > TL
      f_liq_min = (f_liq(1) + f_ice(1) * ro_ice / ro_liq) * f_ell_min;
      f_liq_res = max(f_liq_res, f_liq_min / (1 - f_ice(1)));
   end

   % Compute evaporation/condensation/sublimation.
   [f_ice, f_liq, d_con, x_err] = ICESUBL(f_ice, f_liq, d_con, ...
      ro_ice, ro_liq, Ls, Lv, d_pevp, f_min, f_liq_res);

   % Budget evap / subl.
   d_evp = d_evp + f_liq - xf_liq;

   % Update the top layer temperature, if the node is in the melt zone
   if T(1) > TL % || f_liq(1) > f_liq_min
      T(1) = MELTCURVE(T(1), f_ice(1), f_liq(1), ro_ice, ro_liq, fcp, Tf);
   end

   % Below here:
   % d_XXX = f_liq - xf_liq - d_evp
   % where f_liq on the rhs is f_liq_new
end
