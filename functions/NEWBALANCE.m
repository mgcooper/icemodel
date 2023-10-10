function [Qle,Qh,Qe,Qc,Qm,Qf,balance,Tsfc,Tflag] = NEWBALANCE(Tair,wspd,rh, ...
      Qsi,Qli,albedo,xk_eff,T_old,Tfp,dy_p,chi,xTsfc,xkappa,ro_air,Cp_air,  ...
      emiss_sfc,Stef_Boltz,gravity,xLs,opts,ea,Pa,De_h)
   %NEWBALANCE Compute the surface energy-balance and solve for tsfc

   % Atmospheric vapor pressure from relative humidity data.
   ea = VAPPRESS(rh,Tair,Tfp);

   % Compute the average station pressure.
   Pa = PRESSURE(opts.topo);

   % Compute the turbulent exchange coefficients.
   De_h = EXCOEFS(wspd,xkappa,opts);

   % Compute the flux contribution due to conduction.
   Qc = CONDUCT(xk_eff,T_old,dy_p,xTsfc,Tfp);

   % Solve the energy balance for the surface temperature.
   [Tsfc,Tflag] = SFCTEMP(Tair,Qsi,Qli,ea,albedo,De_h,Pa,wspd,    ...
      ro_air,Cp_air,emiss_sfc,Stef_Boltz,gravity,xLs, ...
      xkappa,Tfp,Qc,chi,opts);

   % Make the Tsfc_0 <= 0 C for surface flux calculations.
   %   Let Tsfc remain > Tf for the upper boundary condition on ICE_ENERGY
   Tsfc0 = MELTTEMP(Tsfc,Tfp);

   % Compute the stability function.
   stability = STABLEFN(Tair,Tsfc0,wspd,gravity,xkappa,opts);

   % Compute the water vapor pressure at the surface.
   es0 = VAPOR(Tsfc0,Tfp);

   % Compute the latent heat flux.
   Qe = LATENT(De_h,stability,ea,es0,ro_air,xLs,Pa);

   % Compute the sensible heat flux.
   Qh = SENSIBLE(De_h,stability,Tair,Tsfc0,ro_air,Cp_air);

   % Compute the longwave flux emitted by the surface.
   Qle = LONGOUT(Tsfc0,emiss_sfc,Stef_Boltz);

   % Compute the energy flux available for melting or freezing.
   [Qm,Qf] = MFENERGY(albedo,Qsi,Qli,Qle,Qh,Qe,Qc,Tsfc,Tfp,  ...
      Tair,wspd,gravity,De_h,ea,ro_air,xLs,Pa,Cp_air, ...
      emiss_sfc,Stef_Boltz,xkappa,xk_eff,T_old,dy_p,  ...
      chi,opts);

   % Perform an energy balance check.
   bal = ENBAL(albedo,Qsi,Qli,Qle,Qh,Qe,Qc,Qm,chi);

   % Compute the effective temperature error
   Terr = TEMP_ERROR(ro_sno,cp_sno,xLs,frac_air,drovdT,bal);

   % For a 'skin' surface energy balance model, reset Tsfc
   if opts.skin_model == true || opts.skin_melt == true
      Tsfc = Tsfc0;
   end
end
