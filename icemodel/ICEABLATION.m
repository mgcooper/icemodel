function [surf_mlt, surf_frz, surf_sub, surf_con, surf_rof] = ICEABLATION( ...
      Qm,Qe,Qf,surf_mlt,surf_frz,surf_rof,surf_sub,surf_con,ro_liq,Lf,Ls,dt, ...
      opts,dz,ro_ice,f_ice,f_liq) %#ok<*INUSD>
   %ICEABLATION Calculate surface melt/freeze, sublimation, and ablation rate
   %
   % note:
   % isubl = -Qe/(Ls*row)*dt [m w.e.]
   % hsubl = isubl*row/roi   [m i.e.]
   % fsubl = hsubl/ht = isubl*row/roi/ht; (see f_ice update in code above)
   %
   % See also:

   % add these back to replicate the original behavior if needed
   % f_ice, f_liq, imelt, isubl

   % All the liquid fluxes are added to surf_runoff. Each individual flux is
   % saved as a cumulative sum. Runoff is calculated as a post-process step
   % (see SURF_RUNOFF), which accounts for condensation and melt available for
   % runoff

   % These are not needed unless they are returned as in the original behavior
   %imelt = 0.0;
   %isubl = 0.0;
   %icond = 0.0;
   %ifreeze = 0.0;
   %runoff = runoff + irain;

   % surface sublimation / deposition
   if Qe < 0.0

      isubl = - Qe / (Ls * ro_liq) * dt; % [m w.e.]
      surf_sub = surf_sub + isubl;

      % NOTE: this happens in ICEMF instead
      %
      % reduce the ice surface by sublimation and melt in the top layer
      %    f_ice(1) = f_ice(1) - isubl*ro_liq/ro_ice/dz(1);
      %
      % previously: ice_surf = ice_surf - isubl * ro_liq / ro_ice;

   elseif Qe > 0.0

      icond = Qe / (Ls * ro_liq) * dt;
      surf_con = surf_con + icond;

      % assume condensation runs off
      %    f_liq(1) = f_liq(1) + icond*ro_liq/ro_ice/dz(1);
      %
      % previously: surf_rof = surf_rof + icond;
   end

   % surface melt / freeze
   if (Qm>0.0 && strcmp(opts.simmodel, 'skinmodel')) % || (Qm>0.0 && opts.skinmelt == true)

      imelt = Qm / (Lf * ro_liq) * dt;
      surf_mlt = surf_mlt + imelt;

      % decrease the ice surface and add imelt to surf_runoff
      surf_rof = surf_rof + imelt;
      % f_ice(1) = f_ice(1) - imelt*ro_liq/ro_ice/dz(1);

      % previously: ice_surf = ice_surf - imelt*ro_liq/ro_ice;

      % this would instead add the melt to the upper grid cell
      % f_liq(1) = f_liq(1) + imelt/dz;

   elseif (Qf > 0.0 && opts.skinfreeze == true)

      ifreeze = Qf / (Lf * ro_liq) * dt;
      surf_frz = surf_frz + ifreeze;
   end
end
