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
   %  icemodel.kernels.buckVaporPressure.
   %
   % See also: icemodel.physicalConstant, VAPORINIT
   %
   %#codegen

   narginchk(0, Inf);

   % -----------------------------------------------------------------------
   % Ambaum (2020) Rankine-Kirchhoff vapor coefficients (canonical)
   % -----------------------------------------------------------------------
   % Computed from physical constants via VAPORINIT. The formula is:
   %   es = a * exp(b / T) * T ^ c   [Pa]

   [al_, bl_, cl_, ai_, bi_, ci_] = VAPORINIT();

   % -----------------------------------------------------------------------
   % Vapor transport parameters
   % -----------------------------------------------------------------------

   params = struct( ...
      'al', al_, ...                  % Rankine-Kirchhoff coefficient over liquid [Pa K^-c]
      'bl', bl_, ...                  % Rankine-Kirchhoff coefficient over liquid [K]
      'cl', cl_, ...                  % Rankine-Kirchhoff exponent over liquid [1]
      'ai', ai_, ...                  % Rankine-Kirchhoff coefficient over ice [Pa K^-c]
      'bi', bi_, ...                  % Rankine-Kirchhoff coefficient over ice [K]
      'ci', ci_, ...                  % Rankine-Kirchhoff exponent over ice [1]
      'nd', 6, ...                    % Vapor diffusivity temperature exponent [1]
      'De0', 9e-5, ...               % Reference vapor diffusivity [m2 s-1]
      'De0_sntherm', 9.2e-5, ...     % SNTHERM89 reference vapor diffusivity [m2 s-1]
      ...
      ... % ----------------------------------------------------------------
      ... % Surface and radiation parameters
      ... % ----------------------------------------------------------------
      'emiss', 0.98, ...             % Surface emissivity for ice [1]
      'S0', 1361.0, ...              % Solar constant [W m-2]
      'n_ice', 1.31, ...            % Refractive index (real) of ice [1]
      'n_liq', 1.33, ...            % Refractive index (real) of water [1]
      'n_air', 1.00, ...            % Refractive index (real) of air [1]
      ...
      ... % ----------------------------------------------------------------
      ... % Empirical and derived atmospheric parameters
      ... % ----------------------------------------------------------------
      'N0', 0.08, ...               % Marshall-Palmer parameter [cm-4]
      'psychro', 66.1, ...          % Psychrometric constant [Pa K-1]
      'dalr', 9.76, ...             % Dry adiabatic lapse rate [K km-1]
      'malr', 5.1, ...              % Moist adiabatic lapse rate [K km-1]
      'scale_ht', 8434.5, ...       % Scale height assuming average temperature [m]
      ...
      ... % ----------------------------------------------------------------
      ... % Model tuning parameters
      ... % ----------------------------------------------------------------
      'fcp', 100, ...               % Freezing curve parameter [K-1]
      ...
      ... % ----------------------------------------------------------------
      ... % Time conventions
      ... % ----------------------------------------------------------------
      'hrsperday', 24, ...          % Hours per day [hr]
      'secperhr', 3600 ...          % Seconds per hour [s]
      );

   % Derived values
   params.emissSB = params.emiss * 5.670374419e-8;  % emissivity x Stefan-Boltzmann [W m-2 K-4]
   params.fcpsq = params.fcp ^ 2;                   % Freezing curve parameter squared [K-2]
   params.secperday = params.hrsperday * params.secperhr; % Seconds per day [s]

   if (nargout == 1 && nargin == 0) || (nargin > 0 && strcmp('all', varargin{1}))
      varargout{1} = params;
      return
   end
   for n = 1:nargin
      arg = varargin{n};
      varargout{n} = params.(arg);
   end
end
