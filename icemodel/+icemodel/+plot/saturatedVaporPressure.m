function saturatedVaporPressure(T)

   if nargin < 1
      T = 230:290;
   end

   % Buck (1981) — archived reference
   [es_1, des_dT_1] = icemodel.kernels.buckVaporPressure(T, 273.16, false);

   % Ambaum (2020) — production (VAPPRESS2 returns full derivative chain)
   [es_2, des_dT_2, d2es_dT2_2, ro_vap_2, dro_vapdT_2, d2ro_vapdT2_2] = ...
      VAPPRESS2(T, false);

   % GETGAMMA comparison
   f_ice = ones(size(T));
   f_liq = zeros(size(T));
   ro_ice = 917;
   [k_liq, Rv, Ls, Tf] = icemodel.physicalConstant("k_liq", "Rv", "Ls", "Tf");
   [k_eff, k_vap] = GETGAMMA(T, f_ice, f_liq, ro_ice, k_liq, Ls, Rv, Tf);

   figure;
   tiledlayout('flow', TileSpacing='compact', Padding='compact')

   nexttile; hold on
   plot(T, es_1);
   plot(T, es_2, ':');
   ylabel('e_s [Pa]')
   set(gca, 'YScale', 'log');
   legend('Buck', 'Ambaum', 'Location', 'northwest')

   nexttile; hold on
   plot(T, des_dT_1);
   plot(T, des_dT_2, ':');
   ylabel('de_s/dT [Pa K^{-1}]')
   set(gca, 'YScale', 'log');
   legend('Buck', 'Ambaum', 'Location', 'northwest')

   nexttile; hold on
   plot(T, d2es_dT2_2);
   ylabel('d^2e_s/dT^2 [Pa K^{-2}]')
   set(gca, 'YScale', 'log');
   legend('Ambaum', 'Location', 'northwest')

   % rho_vap
   nexttile; hold on
   plot(T, ro_vap_2);
   ylabel('\rho_v [kg m^{-3}]')
   set(gca, 'YScale', 'log');
   legend('Ambaum', 'Location', 'northwest')

   % d_rho_vap / dT
   nexttile; hold on
   plot(T, dro_vapdT_2);
   ylabel('d\rho_v / dT [kg m^{-3} K^{-1}]')
   set(gca, 'YScale', 'log');
   legend('Ambaum', 'Location', 'northwest')

   % d2_rho_vap / dT2
   nexttile; hold on
   plot(T, d2ro_vapdT2_2);
   ylabel('d^2\rho_v / dT^2 [kg m^{-3} K^{-2}]')
   set(gca, 'YScale', 'log');
   legend('Ambaum', 'Location', 'northwest')

   figure;
   tiledlayout('flow', TileSpacing='compact', Padding='compact')
   nexttile; hold on
   plot(T, k_eff)
   ylabel('k_{eff} [W m^{-1} K^{-1}]')
   legend('Ambaum', 'Location', 'northwest')

   nexttile; hold on
   plot(T, k_vap)
   ylabel('k_{vap} [W m^{-1} K^{-1}]')
   legend('Ambaum', 'Location', 'northwest')

end
