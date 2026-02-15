
liqflag = false;
tair = 274;
swd = 500;
lwd = 300;
albedo = 0.5;
psfc = 8600;
wspd = 5;
ppt = 0;
tppt = 274;
rh = 80;
ea = VAPPRESS(tair, 273.16, liqflag) .* rh/100;
k_eff = 2;

z_0 = 1e-3;
z_tair = 3;
[De, scoef, wcoef] = WINDCOEF(wspd, z_0, z_tair);

dz = 0.04 * ones(500, 1);

Qc = 0;
chi = 0.5;

Z = cumsum(dz);
T = tair * ones(500, 1);
T = T .* exp(-0.1 * Z./Z(end));

solver_options = [-1, -2, 0, 1, 2];

for solver = solver_options(:)'

   [cv_air, cv_liq, emiss, SB, Tf, roL] = icemodel.physicalConstant(...
      "cv_air", "cv_liq", "emiss", "SB", "Tf", "roLv");
   
   % standard approach (for "slow" code)
   Tsfc = SEBSOLVE(tair, swd, lwd, ...
      albedo, wspd, ppt, tppt, ...
      psfc, De, ea, cv_air, cv_liq, ...
      emiss, SB, Tf, chi, roL, scoef, ...
      liqflag, tair, T, k_eff, ...
      dz, solver);
end
