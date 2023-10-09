function ksnow = GETKTHERMAL(T,f_ice,ro_ice,ro_sno,opts)
%--------------------------------------------------------------------------

% gamma_ice is bulk density of ice i.e. 'dry snow density' [g/cm3]
gam_ice        =   ro_ice.*f_ice;

switch opts.kthermal

   % Sturm et al. 1997 (dry snow from X - X oC, density in g/cm3)
   case 1
      ro       =   gam_ice./1000;
      i        =   ro < 0.156;
      ksnow    =   0.138 - 1.01.*ro + 3.233.*ro.^2;
      ksnow(i) =   0.023 + 0.234.*ro(i);

      % Anderson, 1976 (wet snow from X - X oC, density in g/cm3)
   case 2
      ksnow    =   0.021 + 2.5.*(ro_sno./1000).^2;

      % Aggerwal, (wet snow Himalaya, -3:-14oC, 120:420 kg/m3, density in kg/m3)
   case 3
      a        =   0.00395;
      b        =   0.00084;
      c        =   -1.7756e-6;
      d        =   3.80635e-9;
      ksnow    =   a + b.*gam_ice + c.*gam_ice.^2 + d.*gam_ice.^3;

      % Yen, 1981 (used in HIRHAM), ratio of densities so use consistent units
      % (note: see schwerdtfeger eq. 36 as well)
   case 4
      kice     =   2.22362;
      ksnow    =   kice .* (ro_sno./ro_liq).^1.885;

      % Calonne et al, 2017 (dry snow, density in kg/m3)
   case 5
      a        =   2.5e-6;
      b        =   1.23e-4;
      c        =   0.024;
      ksnow    =   a.*gam_ice.^2 - b.*gam_ice + c;
      ksnow    =   round(ksnow, 3);
      % kcal   =   0.024 - 1.23*10^(-4).*ro_s + 2.5*10^-6 .* ro_s.^2

      % van Dusen, 1929 (lower limit)
   case 6
      a        =   2.1e-2 + 4.2e-4;
      b        =   2.2e-9;
      ksnow    =   a .* ro_sno + b .* ro_sno.^3;

   case 7
      ksnow    =   1.8.*ones(size(ro_sno));

      % Glacier ice in Atarctica
   case 8
      ksnow    =   2.74.*ones(size(ro_sno));

      % Pure ice at 0oC
   case 9
      ksnow    =   2.174.*ones(size(ro_sno));

      % Mixture theory, with temperature dependence based on Schwerdtfeger, 1963
   case 10
      kice     =   9.828 .* exp(-5.7e-3.*T);
      ksnow    =   2 .* kice .* gam_ice ./ (3.*ro_ice./1000-gam_ice./1000);

    % % for low-density would be:
    % s        =   1-1/porosity^3;
    % xk_snow  =   (2+s)*s/((1+s)^2).*k_ice;

      % Calonne etal 2019 Eq.1 & Calonne etal 2011 Eq.12, density in kg/m3
      % NOTE: only valid from 550-917 kg/m3
   case 11
      i        =   gam_ice<550;
      kfirn    =   2.107 + 0.003618.*(gam_ice-ro_ice);
      ksnow    =   kfirn;
      ksnow(i) =   2.107 + 0.003618.*(550-ro_ice);

      % Calonne etal 2019, Eq.5
   case 12
      th       =  1./(1+exp(-0.04.*(gam_ice-450)));
      kfirn    =  2.107 + 0.003618.*(gam_ice-917);               %ref
      ksnow    =  0.024 - 1.23e-4.*gam_ice + 2.5e-6.*gam_ice.^2; %ref
      kiceT    =  9.828 .* exp(-5.7e-3.*T);
      ksnow    =  (1-th).*0.47461.*kiceT.*ksnow+th.*kiceT./2.107.*kfirn;

      % for reference (above is identical but faster):
      %         ki_ref      =   2.107;
      %         ka_ref      =   0.024;
      %         k_airT      =   0.024;
      %         xk_snow     =   (1-th).*(k_iceT.*k_airT)./(ki_ref*ka_ref).*k_snow ...
      %                             + th.*(k_iceT./ki_ref).*k_firn;
      % note: air thermal k temperature dependence is extremely low and for
      % the range of temperatures in near-surface glacier ice is
      % effectively constant 0.024
end