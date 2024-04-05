function varargout = physicalConstant(varargin)
   %PHYSICALCONSTANT return the value of a physical constant.
   %
   % C = physicalConstant(NAME) returns the physical constant C specified by
   % NAME.
   %
   % Multiple constants can be retrieved at once, e.g.,
   %
   % [C1, C2, C3] = physicalConstant('Name1', 'Name2', 'Name3')
   % The returned physical constants use International System of Units (SI).
   %
   % Example:
   %   Obtain the speed of light in the free space and the Boltzmann constant.
   %
   % [c, b] = physicalConstant('LightSpeed', 'Boltzmann')
   %
   %#codegen

   % see: IAPWS release on the Properties of Ice and Supercooled Water.

   narginchk(0,Inf);

   % name = validatestring(name,icemodel.physconlist,'physicalConstant','name', 1);


   % Define physical constants and accepted values
   constants = struct( ...
      'Tf', 273.16, ...             % Freezing point (Triple point), [K]
      'Lv', 2.501e6, ...            % Latent heat of vaporization at 0°C [J kg-1]
      'Lf', 3.3355e5, ...           % Latent heat of fusion at 0°C [J kg-1]
      'Ls', 2.834e6, ...            % Latent heat of sublimation at 0°C [J kg-1]
      'ro_air', 1.293, ...          % Density of air at 0°C and 1 atm [kg m-3]
      'ro_ice', 916.7, ...          % Density of ice at 0°C [kg m-3]
      'ro_liq', 999.8395, ...       % Density of water at 0°C [kg m-3]
      'cp_air', 1005.0, ...         % Specific heat capacity of air at constant pressure at 0°C [J kg-1 K-1]
      'cp_liq', 4218.0, ...         % Specific heat capacity of water at 0°C [J kg-1 K-1]
      'cp_ice', 2093.0, ...         % Specific heat capacity of ice at 0°C [J kg-1 K-1]
      'cp_vap', 1846.0, ...         % specific heat capacity of water vapor at constant pressure [J kg-1 K-1]
      'k_liq', 0.561, ...           % Thermal conductivity of water at 0°C [W m-1 K-1]
      'k_ice', 2.22, ...            % Thermal conductivity of ice at 0°C [W m-1 K-1]
      'Rd', 287.0, ...              % Gas constant for dry air [J kg-1 K-1]
      'Rv', 461.5, ...              % Gas constant for water vapor [J kg-1 K-1]
      'SB', 5.670374419e-8, ...     % Stefan-Boltzmann constant [W m-2 K-4]
      'emiss', 0.98, ...            % Surface emissivity for ice [1]
      'gravity', 9.80665, ...       % Gravitational acceleration [m s-2]
      'kappa', 0.4, ...             % von Karman constant [1]
      'kappa_p', 0.2854, ...        % Poisson constant for dry air [1] (Rd/cp_air)
      'epsilon', 0.622, ...         % Ratio of molar mass of dry air to that of water vapor [1]
      'P0', 101325.0, ...           % One atmosphere [Pa]
      'es0', 611.2, ...             % Reference saturation vapor pressure at 0°C [Pa]
      'S0', 1361.0, ...             % Solar constant [W m-2]
      'N0', 0.08, ...               % Marshall-Palmer parameter [cm-4]
      'psychro', 66.1, ...          % Psychrometric constant [Pa K-1]
      'dalr', 9.76, ...             % Dry adiabatic lapse rate [K km-1]
      'malr', 5.1, ...              % Moist adiabatic lapse rate [K km-1]
      'fcp', 100, ...               % Freezing point depression constant [K-1]
      'scale_ht', 8434.5, ...       % Scale height assuming average temperature [m]
      'hrsperday', 24, ...          % Hours per day [hr]
      'secperhr', 3600 ...          % Seconds per hour [s]
      );

   % Assign derived values based on physical constants

   % Volumetric heat capacities of dry air, ice, and liquid water at constant
   % pressure and 0°C [J m-3 K-1]
   constants.cv_air = constants.ro_air * constants.cp_air;  % [J m-3 K-1]
   constants.cv_ice = constants.ro_ice * constants.cp_ice;  % [J m-3 K-1]
   constants.cv_liq = constants.ro_liq * constants.cp_liq;  % [J m-3 K-1]

   % Volumetric enthalpy of vaporization, sublimation, and fusion
   constants.roLv = constants.ro_air * constants.Lv;        % [J m-3]
   constants.roLs = constants.ro_air * constants.Ls;        % [J m-3]
   constants.roLf = constants.ro_liq * constants.Lf;        % [J m-3]

   % Need to add for VAPORHEAT:
   % Ls / Rv
   % Lv / Rv

   % Water vapor coefficients [kg m-3 K]
   ci = 611.15 * exp(constants.Ls / (constants.Rv * constants.Tf)) / constants.Rv;
   cl = 611.21 * exp(constants.Lv / (constants.Rv * constants.Tf)) / constants.Rv;

   % Ratios of intrinsic ice and water density
   constants.ro_iwe = constants.ro_ice / constants.ro_liq;  % [1]
   constants.ro_wie = constants.ro_liq / constants.ro_ice;  % [1]

   % Freezing point depression parameter ^2
   constants.fcpsq = constants.fcp^2;                       % [K-2]

   % Emissivity x Stefan-Boltzmann constant [W m-2 K-4]
   constants.emissSB = constants.emiss * constants.SB;      % [W m-2 K-4]

   % Time conversions
   constants.secperday = constants.hrsperday * constants.secperhr;

   if (nargout == 1 && nargin == 0) || (strcmp('all', varargin{1}))
      varargout{1} = constants;
   end
   for n = 1:nargin
      arg = varargin{n};
      varargout{n} = constants.(arg);
   end
end

