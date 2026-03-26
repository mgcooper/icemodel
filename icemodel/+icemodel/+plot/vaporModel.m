function vaporModel(T)
   %VAPORMODEL Three-formulation diagnostic comparison plot.
   %
   %  vaporModel() plots saturation vapor pressure and derived quantities for
   %  Buck (1981), Ambaum (2020) / Romps (2021), and the reference
   %  Clausius-Clapeyron kernel over the temperature range 230-290 K.
   %
   %  Figure 1: es, des_dT, d2es_dT2, ro_vap, dro_vapdT, k_vap
   %  Figure 2: Relative differences between Buck and Ambaum
   %
   % See also: VAPPRESS2, icemodel.kernels.buckVaporModel,
   %           icemodel.kernels.saturationVaporPressure

   if nargin < 1
      T = (230:0.5:290)';
   end
   T = T(:);

   Tf = 273.16;
   Ls = icemodel.physicalConstant('Ls');

   % ---------------------------------------------------------------------
   % 1. Buck (1981) — archived reference (all outputs from buckVaporModel)
   % ---------------------------------------------------------------------
   [es_buck, des_dT_buck, ro_vap_buck, dro_vapdT_buck, k_vap_buck] = ...
      icemodel.kernels.buckVaporModel(T, Tf, false);

   % ---------------------------------------------------------------------
   % 2. Ambaum (2020) / Romps (2021) — production
   % ---------------------------------------------------------------------
   [es_amb, des_dT_amb, d2es_dT2_amb, ro_vap_amb, dro_vapdT_amb] = ...
      VAPPRESS2(T, false);

   % Ambaum k_vap
   [De0, nd] = icemodel.parameterLookup('De0', 'nd');
   De_amb = De0 * (T / Tf) .^ nd;
   k_vap_amb = Ls * De_amb .* dro_vapdT_amb;

   % ---------------------------------------------------------------------
   % 3. Reference Clausius-Clapeyron kernel (Romps exact, es only)
   % ---------------------------------------------------------------------
   es_cc = icemodel.kernels.saturationVaporPressure(T, false);

   % =====================================================================
   % Figure 1: Absolute quantities
   % =====================================================================
   figure('Name', 'Vapor Model: Three Formulations', ...
      'Position', [100 100 1200 800]);
   tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
   title(tl, 'Saturation vapor quantities (over ice)', 'FontWeight', 'bold')

   % -- es --
   nexttile; hold on
   plot(T, es_buck, 'b-', 'LineWidth', 1.2)
   plot(T, es_amb, 'r--', 'LineWidth', 1.2)
   plot(T, es_cc, 'k:', 'LineWidth', 1.2)
   set(gca, 'YScale', 'log')
   ylabel('e_s [Pa]')
   xlabel('T [K]')
   legend('Buck', 'Ambaum', 'C-C ref', 'Location', 'northwest')
   title('Saturation vapor pressure')

   % -- des_dT --
   nexttile; hold on
   plot(T, des_dT_buck, 'b-', 'LineWidth', 1.2)
   plot(T, des_dT_amb, 'r--', 'LineWidth', 1.2)
   set(gca, 'YScale', 'log')
   ylabel('de_s/dT [Pa K^{-1}]')
   xlabel('T [K]')
   legend('Buck', 'Ambaum', 'Location', 'northwest')
   title('First derivative of e_s')

   % -- d2es_dT2 --
   nexttile; hold on
   plot(T, d2es_dT2_amb, 'r--', 'LineWidth', 1.2)
   % Buck d2es_dT2 from es and des_dT (not in buckVaporModel)
   bi_buck = 22.452; ci_buck = 272.55;
   Td = T - Tf;
   d2es_dT2_buck = es_buck .* bi_buck .* ci_buck .* ...
      (bi_buck .* ci_buck - 2 * (ci_buck + Td)) ./ (ci_buck + Td).^4;
   plot(T, d2es_dT2_buck, 'b-', 'LineWidth', 1.2)
   set(gca, 'YScale', 'log')
   ylabel('d^2e_s/dT^2 [Pa K^{-2}]')
   xlabel('T [K]')
   legend('Ambaum', 'Buck', 'Location', 'northwest')
   title('Second derivative of e_s')

   % -- ro_vap --
   nexttile; hold on
   plot(T, ro_vap_buck, 'b-', 'LineWidth', 1.2)
   plot(T, ro_vap_amb, 'r--', 'LineWidth', 1.2)
   set(gca, 'YScale', 'log')
   ylabel('\rho_v [kg m^{-3}]')
   xlabel('T [K]')
   legend('Buck', 'Ambaum', 'Location', 'northwest')
   title('Saturation vapor density')

   % -- dro_vapdT --
   nexttile; hold on
   plot(T, dro_vapdT_buck, 'b-', 'LineWidth', 1.2)
   plot(T, dro_vapdT_amb, 'r--', 'LineWidth', 1.2)
   set(gca, 'YScale', 'log')
   ylabel('d\rho_v/dT [kg m^{-3} K^{-1}]')
   xlabel('T [K]')
   legend('Buck', 'Ambaum', 'Location', 'northwest')
   title('Vapor density derivative')

   % -- k_vap --
   nexttile; hold on
   plot(T, k_vap_buck, 'b-', 'LineWidth', 1.2)
   plot(T, k_vap_amb, 'r--', 'LineWidth', 1.2)
   set(gca, 'YScale', 'log')
   ylabel('k_{vap} [W m^{-1} K^{-1}]')
   xlabel('T [K]')
   legend('Buck', 'Ambaum', 'Location', 'northwest')
   title('Vapor thermal diffusion coeff.')

   % =====================================================================
   % Figure 2: Relative differences (Buck vs Ambaum)
   % =====================================================================
   figure('Name', 'Vapor Model: Relative Differences', ...
      'Position', [150 80 1200 500]);
   tl2 = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
   title(tl2, 'Relative difference: (Buck - Ambaum) / Ambaum', ...
      'FontWeight', 'bold')

   % -- es relative difference --
   nexttile; hold on
   plot(T, (es_buck - es_amb) ./ es_amb * 100, 'b-', 'LineWidth', 1.2)
   plot(T, (es_cc - es_amb) ./ es_amb * 100, 'k:', 'LineWidth', 1.2)
   ylabel('Relative difference [%]')
   xlabel('T [K]')
   yline(0, '--', 'Color', [0.5 0.5 0.5])
   legend('Buck vs Ambaum', 'C-C ref vs Ambaum', 'Location', 'best')
   title('e_s')

   % -- des_dT relative difference --
   nexttile; hold on
   plot(T, (des_dT_buck - des_dT_amb) ./ des_dT_amb * 100, ...
      'b-', 'LineWidth', 1.2)
   ylabel('Relative difference [%]')
   xlabel('T [K]')
   yline(0, '--', 'Color', [0.5 0.5 0.5])
   title('de_s/dT')

   % -- k_vap relative difference --
   nexttile; hold on
   plot(T, (k_vap_buck - k_vap_amb) ./ k_vap_amb * 100, ...
      'b-', 'LineWidth', 1.2)
   ylabel('Relative difference [%]')
   xlabel('T [K]')
   yline(0, '--', 'Color', [0.5 0.5 0.5])
   title('k_{vap}')
end
