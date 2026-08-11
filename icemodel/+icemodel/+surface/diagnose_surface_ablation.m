function [surf_mlt, surf_frz, surf_sub, surf_con, surf_rof] = diagnose_surface_ablation( ...
      Qm, Qe, Qf, surf_mlt, surf_frz, surf_rof, surf_sub, surf_con, dt, opts)
   %DIAGNOSE_SURFACE_ABLATION Diagnose cumulative surface ablation terms.
   %
   % The sublimation terms are:
   % isubl = -Qe/(Ls*row)*dt [m w.e.]
   % hsubl = isubl*row/roi   [m i.e.]
   % fsubl = hsubl/ht = isubl*row/roi/ht; (see f_ice update in code above)
   % The inputs are the cumulative state variables and the fluxes of the
   % current timestep. This function reads its own physical constants, so the
   % caller passes only dt and opts.
   %
   % See also:
   %
   %#ok<*INUSD>
   %#codegen

   persistent ro_liq Lf Ls
   if isempty(ro_liq)
      [ro_liq, Lf, Ls] = icemodel.physicalConstant('ro_liq', 'Lf', 'Ls');
   end

   % add these back to replicate the original behavior if needed
   % f_ice, f_liq, imelt, isubl

   % This function adds every liquid flux to surf_runoff and saves each flux
   % as a cumulative sum. A post-processing step computes runoff (see
   % SURF_RUNOFF). That step accounts for the condensation and the melt that
   % are available for runoff.

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

      % NOTE: this happens in budget_surface_mass_balance instead
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
   if (Qm>0.0 && strcmp(opts.smbmodel, 'skinmodel')) % || (Qm>0.0 && opts.skinmelt == true)

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
