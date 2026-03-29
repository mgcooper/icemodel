function varargout = parameterLookup(varargin)
   %PARAMETERLOOKUP Return the value of a model parameter.
   %
   %  P = parameterLookup(NAME) returns the model parameter P specified by
   %  NAME.
   %
   %  Multiple parameters can be retrieved at once, e.g.,
   %
   %  [P1, P2, P3] = parameterLookup('Name1', 'Name2', 'Name3')
   %
   %  This function is the canonical source of tunable model parameters,
   %  empirical coefficients, and formulation-specific derived quantities.
   %  True physical constants belong in icemodel.physicalConstant.
   %
   %  Canonical vapor coefficients are from Ambaum (2020), computed by
   %  VAPORINIT. Buck (1981) coefficients are preserved in
   %  icemodel.kernels.buckVaporModel.
   %
   % See also: icemodel.physicalConstant, VAPORINIT
   %
   %#codegen

   narginchk(0, Inf);

   % -----------------------------------------------------------------------
   % Ambaum (2020) / Romps (2021) Rankine-Kirchhoff vapor coefficients
   % -----------------------------------------------------------------------
   % Computed from physical constants via VAPORINIT. The formula is:
   %   es = a * exp(b / T) * T ^ c   [Pa]

   [al_, bl_, cl_, ai_, bi_, ci_] = VAPORINIT();

   % Define the parameter struct.
   params = struct( ...
      ... ----------------------------------------------------------------
      ... Vapor transport parameters
      ... ----------------------------------------------------------------
      'al', al_, ...    Rankine-Kirchhoff coefficient over liquid [Pa K^-c]
      'bl', bl_, ...    Rankine-Kirchhoff coefficient over liquid [K]
      'cl', cl_, ...    Rankine-Kirchhoff exponent over liquid [1]
      'ai', ai_, ...    Rankine-Kirchhoff coefficient over ice [Pa K^-c]
      'bi', bi_, ...    Rankine-Kirchhoff coefficient over ice [K]
      'ci', ci_, ...    Rankine-Kirchhoff exponent over ice [1]
      'De0', 9e-5, ...  Reference vapor diffusivity at 1000mb and 273.15k [m2 s-1]
      'nd', 6, ...      Vapor diffusivity temperature exponent [1]
      ...
      ... f_liq threshold above which thermodynamics are wrt liquid phase
      'f_liq_phase_switch_threshold', 0.02, ...
      ...
      ... ----------------------------------------------------------------
      ... Grain growth parameters (Jordan 1991, SNTHERM89 Eqs. 33-34)
      ... ----------------------------------------------------------------
      'g1', 5.0e-7, ...       Dry snow grain growth constant [m4 kg-1]
      'g2', 4.0e-12, ...      Wet snow grain growth constant [m2 s-1]
      'r_max', 2.5e-3, ...    Maximum grain radius [m]
      'Uv_max', 1.0e-6, ...   Max vapor flux for grain growth [kg m-2 s-1]
      ...
      ... ----------------------------------------------------------------
      ... Calonne (2019) thermal conductivity (THERMALK)
      ... ----------------------------------------------------------------
      'cal_ki_ref', 2.107, ...      Reference ice conductivity at -3°C [W m-1 K-1]
      'cal_ka_ref', 0.024, ...      Air thermal conductivity at -3°C [W m-1 K-1]
      'cal_rho_ref', 917.0, ...     Reference ice density [kg m-3]
      'cal_rho_mid', 450.0, ...     Snow/firn logistic midpoint density [kg m-3]
      'cal_a_theta', 0.04, ...      Logistic transition steepness [m3 kg-1]
      'cal_c_firn', 0.003618, ...   Firn density slope [W kg-1 m-2]
      'cal_c_sno1', -0.000123, ...  Snow linear density coefficient [W kg-1 m-2]
      'cal_c_sno2', 2.5e-6, ...     Snow quadratic density coefficient [W kg-2 m-1]
      'cal_c_kT', 9.828, ...        Temperature-dependent ice k coefficient [W m-1 K-1]
      'cal_c_kT_exp', 0.0057, ...   Temperature-dependent ice k exponent [K-1]
      ...
      ... ----------------------------------------------------------------
      ... Surface parameters
      ... ----------------------------------------------------------------
      'emiss', 0.98, ...            Surface emissivity for ice [1]
      ...
      ... ----------------------------------------------------------------
      ... Empirical and derived atmospheric parameters
      ... ----------------------------------------------------------------
      'N0', 0.08, ...               Marshall-Palmer parameter [cm-4]
      'psychro', 66.1, ...          Psychrometric constant [Pa K-1]
      'dalr', 9.76, ...             Dry adiabatic lapse rate [K km-1]
      'malr', 5.1, ...              Moist adiabatic lapse rate [K km-1]
      'scale_ht', 8434.5, ...       Scale height assuming average temperature [m]
      ...
      ... ----------------------------------------------------------------
      ... Model tuning parameters
      ... ----------------------------------------------------------------
      'fcp', 100, ...               Freezing curve slope parameter [K-1]
      ...                           100 approximates a step function
      ... ----------------------------------------------------------------
      ... Time conventions
      ... ----------------------------------------------------------------
      'hrsperday', 24, ...          Hours per day [hr]
      'secperhr', 3600 ...          Seconds per hour [s]
      );

   % Derived values
   SB = icemodel.physicalConstant('SB');
   params.emissSB = params.emiss * SB;
   params.fcpsq = params.fcp ^ 2;
   params.secperday = params.hrsperday * params.secperhr;

   % Parse outputs
   if (nargout == 1 && nargin == 0) || (nargin > 0 && strcmp('all', varargin{1}))
      varargout{1} = params;
      return
   end
   for n = 1:nargin
      arg = varargin{n};
      varargout{n} = params.(arg);
   end
end
