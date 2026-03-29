function k_sno = thermal_conductivity_snow(T, f_ice, ro_ice, ro_sno, reference)
   %THERMAL_CONDUCTIVITY_SNOW Archived multi-option snow conductivity helper.
   %
   %  k_sno = thermal_conductivity_snow(T, f_ice, ro_ice, ro_sno, reference)
   %
   % Description:
   %  Returns snow or porous-ice thermal conductivity k_sno [W m-1 K-1]
   %  using one of the archived literature/reference selectors listed
   %  below. This is the archived multi-option thermal conductivity
   %  function. Production code uses THERMALK (Calonne 2019 Eq. 5 only).
   %
   % Inputs:
   %  T      - Temperature [K]
   %  f_ice  - Ice volume fraction [1]
   %  ro_ice - Intrinsic density of ice [kg m-3]
   %  ro_sno - Bulk snow density [kg m-3]
   %  reference - Conductivity reference selector
   %
   % Options:
   %  "sturm_1997"          - Sturm et al. (1997) dry-snow density fit
   %  "anderson_1976"       - Anderson (1976) wet-snow density fit
   %  "aggarwal_2009"       - Aggarwal et al. (2009) wet-snow density fit
   %  "yen_1981"            - Yen (1981) density-ratio fit
   %  "calonne_2017"        - Calonne et al. (2017) dry-snow density fit
   %  "van_dusen_1929"      - van Dusen (1929) lower-limit density fit
   %  "schwerdtfeger_1963"  - Schwerdtfeger (1963) mixture-theory form
   %                          with temperature-dependent ice conductivity
   %  "calonne_2019_eq1"    - Calonne et al. (2019) high-density firn
   %                          branch, with the low-density clamp retained
   %  "calonne_2019_eq5"    - Calonne et al. (2019) Eq. 5 logistic blend
   %
   % Notes:
   %  g_ice = ro_ice .* f_ice is the dry snow density used by several of the
   %  archived parameterizations below.
   %
   % See also: THERMALK, BULKTHERMALK, thermal_conductivity_firn
   %
   %#codegen

   arguments
      T {mustBeNumeric}
      f_ice {mustBeNumeric}
      ro_ice {mustBeNumeric}
      ro_sno {mustBeNumeric}
      reference (1, 1) string {mustBeMember(reference, [ ...
         "sturm_1997", "anderson_1976", "aggarwal_2009", "yen_1981", ...
         "calonne_2017", "van_dusen_1929", "schwerdtfeger_1963", ...
         "calonne_2019_eq1", "calonne_2019_eq5"])} = "calonne_2019_eq5"
   end

   persistent ro_liq
   if isempty(ro_liq)
      ro_liq = icemodel.physicalConstant('ro_liq');
   end

   % g_ice is bulk density of ice i.e. 'dry snow density' [g/cm3]
   g_ice = ro_ice .* f_ice;

   switch reference

      case "sturm_1997"
         % Sturm et al. 1997
         % Dry snow from X - X oC, density in g/cm3
         ro = g_ice ./ 1000;
         i = ro < 0.156;
         k_sno = 0.138 - 1.01 .* ro + 3.233 .* ro .^ 2;
         k_sno(i) = 0.023 + 0.234 .* ro(i);

      case "anderson_1976"
         % Anderson, 1976
         % Wet snow from X - X oC, density in g/cm3
         k_sno = 0.021 + 2.5 .* (ro_sno ./ 1000) .^ 2;

      case "aggarwal_2009"
         % Aggarwal et al., 2009
         % Wet snow Himalaya, -3:-14oC, 120:420 kg/m3, density in kg/m3
         a = 0.00395;
         b = 0.00084;
         c = -1.7756e-6;
         d = 3.80635e-9;
         k_sno = a + b .* g_ice + c .* g_ice .^ 2 + d .* g_ice .^ 3;

      case "yen_1981"
         % Yen, 1981 (used in HIRHAM), ratio of densities so use consistent
         % units (note: see schwerdtfeger eq. 36 as well)
         kice = 2.22362;
         k_sno = kice .* (ro_sno ./ ro_liq) .^ 1.885;

      case "calonne_2017"
         % Calonne et al, 2017
         % Dry snow, density in kg/m3
         a = 2.5e-6;
         b = 1.23e-4;
         c = 0.024;
         k_sno = a .* g_ice .^ 2 - b .* g_ice + c;
         k_sno = round(k_sno, 3);
         % kcal = 0.024 - 1.23*10^(-4).*ro_s + 2.5*10^-6 .* ro_s.^2

      case "van_dusen_1929"
         % van Dusen, 1929 (lower limit)
         a = 2.1e-2;
         b = 4.2e-4;
         c = 2.2e-9;
         k_sno = a + b .* ro_sno + c .* ro_sno .^ 3;

      case "schwerdtfeger_1963"
         % Schwerdtfeger, 1963, mixture theory w/temperature dependence
         kice = 9.828 * exp(-0.0057 * T); % Yen, 1981
         k_sno = 2 .* kice .* g_ice ./ (3 .* ro_ice - g_ice);

         % % for low-density would be:
         % s = 1-1/porosity^3;
         % k_sno = (2+s)*s/((1+s)^2).*k_ice;

      case "calonne_2019_eq1"
         % Calonne etal 2019 Eq.1 & Calonne etal 2011 Eq.12, density in kg/m3
         % NOTE: only valid from 550-917 kg/m3
         i = g_ice < 550;
         kfirn = 2.107 + 0.003618 .* (g_ice - ro_ice);
         k_sno = kfirn;
         k_sno(i) = 2.107 + 0.003618 .* (550 - ro_ice);

      case "calonne_2019_eq5"
         % Calonne etal 2019, Eq.5
         th = 1 ./ (1 + exp(-0.04 * (g_ice - 450)));
         kfirn = 2.107 + 0.003618 * (g_ice - 917);                % kfirn_ref
         k_sno = 0.024 - 1.23e-4 * g_ice + 2.5e-6 * g_ice .^ 2;   % ksnow_ref
         kiceT = 9.828 * exp(-0.0057 * T);
         k_sno = (1 - th) .* 0.47461 .* kiceT .* k_sno ...
            + th .* kiceT ./ 2.107 .* kfirn;

         % for reference (above is identical but faster):
         % ki_ref = 2.107;
         % ka_ref = 0.024;
         % k_airT = 0.024;
         % xk_snow = (1-th).*(k_iceT.*k_airT)./(ki_ref*ka_ref).*k_snow ...
         %    + th.*(k_iceT./ki_ref).*k_firn;
         % note: air thermal k temperature dependence is extremely low and for
         % the range of temperatures in near-surface glacier ice is
         % effectively constant 0.024

      otherwise
         error('Unsupported snow conductivity reference: %s', reference);
   end
end
