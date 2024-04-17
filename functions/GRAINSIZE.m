function diam = GRAINSIZE(JJ, dt, ro_sno, ...
      f_liq, T, dz, Tf, press, F_lay, ro_ice)

   % This code follows Jordan 1991 (Tech Doc for SNTHERM.89; this is
   %   referenced as "SNTH89" in the comments below).  The equations
   %   that are solved are those in the original SNTHERM89 code.
   %
   % De0          -  Effective diffusion coefficient for water vapor in
   %                 snow (m^2/s) at 100000 Pa and 273.15K
   % press        -  Barometric pressure (Pa)
   % dc(JJ)       -  Vapor diffusion coefficient for = de*change in
   %                 saturation vapor pressure with temperature
   % F_vap(JJ)    -  Average diffusive flux at nodal centroid
   %                 [(upper+lower)/2] (kg/m^2 s)
   % T(JJ)        -  Temperature in center (node) of layer (K)
   % g1           -  Grain growth constant for dry snow (5.0e-7)
   % g2           -  Grain growth constant for wet snow (4.0e-12)
   % diam         -  Nodal effective grain diameter (m)
   % f_liq        -  Volume fraction of liquid water
   % dt           -  Time step (s)
   % e0           -  Saturation vapor pressure at 0 degrees C (mb)
   % Rw           -  Gas constant for water vapor (461.5) (j/kg-K)
   % ro_ice       -  Ice density (917.0) (kg/m^3)
   % ro_sno(JJ)   -  Nodal bulk snow density (kg/m^3)

   % These are the values at the control volume boundaries.
   %  T_layer_bndry(nz_max+1)
   %  dz_bndry(nz_max+1)
   %  diff_coef_bndry(nz_max+1)

   g1 = 5.0e-7;
   g2 = 4.0e-12;
   De0 = 9.2e-5;

   % Calculate the vapor diffusion constants for water(w) and ice(i).
   % e0 /613.6/
   % Rw /461.5/
   % xLvw /2.505e6/
   % xLvi /2.838e6/
   % xLvw_o_Rw = xLvw / Rw
   % c1w = e0 * exp(xLvw_o_Rw / Tf) / Rw
   % xLvi_o_Rw = xLvi / Rw
   % c1i = e0 * exp(xLvi_o_Rw / Tf) / Rw

   Lvw_o_Rw = 5427.952;
   Lvi_o_Rw = 6149.513;

   c1w = 5.6739e8;
   c1i = 7.9639e9;

   pi = 2.0 * acos(0.0);

   % This 'grain_scale' controls allows the user additional control on the
   % grain-growth calculations. Values greater than 1.0 grows grains faster,
   % and values less than 1.0 grows grains slower.
   % grain_scale = 1.0;
   % grain_scale = 0.5;
   grain_scale = 2.5;

   % Initialize the flux output variable so the no-layer outputs make sense.
   for k = 1:nz_max
      F_lay(k) = nan;
   end

   phi = nan(JJ, 1);
   diam = nan(JJ, 1);
   C_kT = nan(JJ, 1);
   D_vap = nan(JJ, 1);
   F_vap = nan(JJ, 1); % vapor flux at nodes
   F_lay = nan(JJ, 1); % vapor flux at upper layer interfaces
   
   
   % Define the snow porosity.
   for k = 1:JJ
      phi(k) = 1.0 - (ro_sno(k) / ro_ice);

      % The only place this value should have problems is if ro_layer is
      % undef. It should be well-constrained before getting to this subroutine.
      if phi(k) < 0.0
         phi(k) = 0.0;
      end
   end

   % Diffusion coefficient for each snow layer.
   for k = 1:JJ

      % Eqn. 20 of SNTH89.
      if f_liq(k) > 0.02

         C_kT(k) = c1w * exp(- Lvw_o_Rw / T(k)) ...
            * (Lvw_o_Rw / T(k) - 1.0) / T(k) ^ 2;
      else

         C_kT(k) = c1i * exp(- Lvi_o_Rw / T(k)) ...
            * (Lvi_o_Rw / T(k) - 1.0) / T(k) ^ 2;
      end

      % SNTH89 assumed the porosity does not affect the fluxes.  SNTH89 had a
      % note that said "The diffusive vapor flux in snow is customarily taken as
      % independent of porosity, which is generally a consequence of the
      % 'hand-to-hand' process of vapor diffusion."  Conceptually, it seems like
      % if the porosity goes to zero the fluxes should stop.  I do that here.
      % If you don't like this idea, you can just comment out these 3 lines.
      % Scale the flux by the available vapor.
      vaporvol = phi(k) - f_liq(k);
      vaporvol = max(0.0, vaporvol);
      C_kT(k) = vaporvol * C_kT(k);

      % Left part of Eqn. 21 of SNTH89.
      Des = De0 * (100000.0 / press) * (T(k) / Tf) ^ 6;

      % Eqn. 21 of SNTH89. The vapor diffusion coefficents.
      D_vap(k) = Des * C_kT(k);

      if D_vap(k) <= 0.0
         % shouldn't happen
      end

   end

   % Include the boundary information in the required arrays.  This information,
   % and the flux calcuations below, follow Glen's 3D thermal sea ice growth
   % model code (see /nice/heat/temp_3d.f, or SeaIce-3D).

   % Number of shifted layers.
   nzm = JJ + 1;

   % Control volume size.
   dz_pm = [0; dz; 0];

   % Temperature.
   T_pm = [T(1); T; T(JJ)];

   % Diffusion coefficients.
   D_pm = [D_vap(1); D_vap; D_vap(end)];

   % Calculate the vapor flux across the control volume walls (Eqn. 21 of
   % SNTH89).  This is: flux = - diff_coef * dT/dz, where diff_coef = Des *
   % C_kT.  See page 45 of Patankar (1980) for a description of what is being
   % done with the harmonic mean calculations.  Here I am solving Patankar's
   % Eqn. 4.8, with delx_e- = 0.5*dx_P, and delx_e+ = 0.5*dx_E.

   for k = 2:nzm

      if dz_pm(k-1) == 0.0 && dz_pm(k) == 0.0

         Uv_bot = 0.0;

      elseif dz_pm(k) == 0.0 && dz_pm(k+1) == 0.0

         Uv_top = 0.0;
      else

         Uv_bot = - 2.0 * (T_pm(k) - T_pm(k-1)) ...
            / (dz_pm(k-1) / D_pm(k-1) ...
            + dz_pm(k) / D_pm(k));

         Uv_top = - 2.0 * (T_pm(k+1) - T_pm(k)) ...
            / (dz_pm(k) / D_pm(k) ...
            + dz_pm(k+1) / D_pm(k+1));
      end

      % Assume the flux at the center of the control volume is the average of
      % the fluxes at the control volume walls.  Also note that we don't care
      % about the direction of the fluxes; we are just assuming that the vapor
      % transport is growing grains, regardless of the direction.  This is used
      % in the grain growth algorithm.
      F_vap(k-1) = (abs(Uv_bot) + abs(Uv_top)) / 2.0;

      % Save a record of the layer fluxes, including the direction of the flow.
      % Here I am saving the flux across the top of each layer. Under the
      % assumption that the flux across the bottom of the bottom layer is zero,
      % this should be enough information to calculate d_flux/d_z and get the
      % mass loss and/or gain in each layer.
      F_lay(k-1) = Uv_top;
   end

   % Because the zero temperature gradient assumed at the boundaries is not
   % realistic, set the boundary fluxes equal to those just inside the
   % boundaries.
   F_vap(1) = F_vap(2);
   F_vap(JJ) = F_vap(JJ-1);

   % Because below we don't allow the abs(fluxes) to be over 1.0e-6,
   % do the same thing here.
   F_lay(F_lay < -1.0e-6) = -1.0e-6;
   F_lay(F_lay > 1.0e-6) = 1.0e-6;
   
   % Convert these to fluxes per dt (instead of per sec). Without this the
   % values are something like 10^-11. 
   F_lay = dt * F_lay;

   % Update the snow grain diameter.
   for k = 1:JJ

      if diam(k) <= 0.0

         error(['execution halted because diam_layer <= 0.0, ' ...
            'layer = %s, diam_layer(k) = %s'],k , diam(k))
      end

      % Dry snow: The cut-off between dry and wet snow is arbitrarily 1e-4

      if f_liq(k) < 1.0e-4

         % The max vapor flux available for growth is arbitrarily set at 1.0e-6.
         % This can be increased if you want to allow larger grains to grow, if
         % the temperature gradients are available to drive greater fluxes.
         % Here I have also included a scaling term that can be used to increase
         % or decrease the growth rate of the grains to better match any
         % observational datasets you might have.  Eqn. 33 of SNTH89.

         if abs(F_vap(k)) < 1.0e-6
            diam(k) = diam(k) + grain_scale * dt * g1 ...
               * abs(F_vap(k)) / diam(k);
         else
            diam(k) = diam(k) + grain_scale * dt * g1 ...
               * 1.0e-6 / diam(k);
         end

         % Wet snow: Different equations for liquid volume fraction
         %   above and below 0.09.
      else
         if f_liq(k) < 0.09
            % Eqn. 34a of SNTH89.
            diam(k) = diam(k) + grain_scale * dt * g2 ...
               * (f_liq(k) + 0.05) / diam(k);-
         else
            % Eqn. 34b of SNTH89.
            diam(k) = diam(k) + grain_scale * dt * g2 ...
               * 0.14 / diam(k);
         end
      end

      % Max grain size set at 5mm, based on Arctic Alaska depth hoar
      %   observations (larger sizes are possible but not common).
      if diam(k) > 5.0e-3
         diam(k) = 5.0e-3;
      end
   end
end
