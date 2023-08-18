function [ro_sno, cp_sno, liqflag, roL, xT, xTsfc, xf_liq, xf_ice, dt_sum, ...
   dt_new, dt_flag] = UPDATESUBSTEP(f_ice, f_liq, ro_ice, ro_liq, ro_air, ...
   cv_ice, cv_liq, T, Tsfc, dt, dt_sum, dt_new, roLv, roLs, dt_min, TINY)

% UPDATE DENSITY, HEAT CAPACITY, DIFFUSION LENGTH SCALE
ro_sno = f_ice.*ro_ice + f_liq.*ro_liq + (1.0 - f_liq - f_ice).*ro_air;
cp_sno = (cv_ice.*f_ice + cv_liq.*f_liq) ./ ro_sno;

% top node contains >2% liquid water
liqflag = f_liq(1) > 0.02;
if liqflag
   roL = roLv;  % ro_air*Lv
else
   roL = roLs;  % ro_air*Ls
end

% STORE PAST VALUES
xT = T;
xTsfc = Tsfc;
xf_liq = f_liq;
xf_ice = f_ice;

% ALLOCATE THIS SUBSTEP TO THE TIMESTEP
dt_sum = dt_sum + dt_new;

% first is true if step incomplete, second if overage will occur
dt_flag = (dt - dt_sum) > TINY && (dt_sum + dt_new - dt) > TINY;

if dt_flag
   dt_new = max(dt - dt_sum, dt_min);
end

% zD = sqrt(k_eff(1)*dt/(ro_sno(1)*cp_sno(1)));
% if zD > dz(1)
%  % placeholder
% end

% an older option used to interpolate between timesteps
% itime = itime + seconds(dt_new);
% this was initialized in INITTIMESTEPS
% itime = Time(1); 