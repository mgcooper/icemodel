function [T, f_liq, f_ice, k_eff, ro_sno, cp_sno] = UPDATESUBL(T, f_liq, ...
      f_ice, ro_ice, ro_liq, ro_air, ro_iwe, cp_ice, cp_liq, k_liq, Qe, ...
      Lv, Ls, Rv, Tf, fcp, dt_new, dz_therm, f_min, liqflag)

   % these have dimensions of f_liq and f_ice i.e. h_liq/h and h_ice/h
   if liqflag == true
      f_liq(1) = min(max(f_liq(1)+Qe/(Lv*ro_liq)*dt_new/dz_therm,0),1);
   else
      f_ice(1) = min(max(f_ice(1)+Qe/(Ls*ro_ice)*dt_new/dz_therm,f_min),1);
   end

   % update the upper temperature
   f_wat = f_liq(1)+f_ice(1).*ro_iwe;
   T(1) = Tf-sqrt((f_wat./f_liq(1)-1))./fcp;

   % update coefficients
   [k_eff, ro_sno, cp_sno] = UPDATESTATE(T,f_ice,f_liq,ro_ice,ro_liq, ...
      ro_air,cp_ice,cp_liq,k_liq,Ls,Rv,Tf);
end
