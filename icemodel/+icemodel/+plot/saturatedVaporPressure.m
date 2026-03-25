function saturatedVaporPressure(T)

   if nargin < 1
      T = 230:290;
   end

   figure

   [es_1, des_dT_1, d2es_dT2_1, ro_vap_1, dro_vapdT_1, d2ro_vapdT2_1] = ...
      VAPPRESS(T, 273.16, true);
   [es_2, des_dT_2, d2es_dT2_2, ro_vap_2, dro_vapdT_2, d2ro_vapdT2_2] = ...
      VAPPRESS2(T, true);

   % Can I get
   f_ice = 1.0;
   f_liq = 0.0;
   ro_ice = 917;
   [k_liq, Rv, Ls, Tf] = icemodel.physicalConstant("k_liq", "Rv", "Ls", "Tf");

   [k_eff, k_vap, dk_eff_dT] = GETGAMMA(T, f_ice, f_liq, ro_ice, k_liq, Ls, Rv, Tf);

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
   plot(T, d2es_dT2_1);
   plot(T, d2es_dT2_2, ':');
   ylabel('d^2e_s/dT^2 [Pa K^{-2}]')
   set(gca, 'YScale', 'log');
   legend('Buck', 'Ambaum', 'Location', 'northwest')

   % rho_vap
   nexttile; hold on
   plot(T, ro_vap_1);
   plot(T, ro_vap_2, ':');
   ylabel('\rho_v [kg m^{-3}]')
   set(gca, 'YScale', 'log');
   legend('Buck', 'Ambaum', 'Location', 'northwest')

   % d_rho_vap / dT
   nexttile; hold on
   plot(T, dro_vapdT_1);
   plot(T, dro_vapdT_2, ':');
   ylabel('d\rho_v / dT [kg m^{-3} K^{-1}]')
   set(gca, 'YScale', 'log');
   legend('Buck', 'Ambaum', 'Location', 'northwest')

   % d2_rho_vap / dT2
   nexttile; hold on
   plot(T, d2ro_vapdT2_1);
   plot(T, d2ro_vapdT2_2, ':');
   ylabel('d^2\rho_v / dT^2 [kg m^{-3} K^{-2}]')
   set(gca, 'YScale', 'log');
   legend('Buck', 'Ambaum', 'Location', 'northwest')


   figure;
   plot(T, dk_eff_dT)
   ylabel('dk_v / dT [W m^{-1} K^{-2}]')
   %    set(gca, 'YScale', 'log');
   legend('Buck', 'Location', 'northwest')
   %    legend('Buck', 'Ambaum', 'Location', 'northwest')

end
