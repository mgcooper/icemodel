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

   narginchk(0,Inf);

   % name = validatestring(name,icemodel.physconlist,'physicalConstant','name', 1);

   % Define all constants
   constants = struct( ...
      'Lv', 2.500e6, ...         % latent heat vapor,              [J kg-1]
      'Lf', 3.34e5, ...          % latent heat fusion,             [J kg-1]
      'Ls', 2.83e6, ...          % latent heat subl.               [J kg-1]
      'Tf', 273.15, ...          % freezing point,                 [K]
      'ro_air', 1.275, ...       % density air,                    [kg m-3]
      'ro_ice', 917.0, ...       % density pure ice,               [kg m-3]
      'ro_liq', 1000.0, ...      % density water,                  [kg m-3]
      'cp_air', 1004, ...        % sp. heat cap. air,              [J kg-1 K-1]
      'cp_liq', 4180.0, ...      % sp. heat cap. water,            [J kg-1 K-1]
      'cp_ice', 2102.0, ...      % sp. heat cap. ice,              [J kg-1 K-1]
      'k_liq', 0.552, ...        % therm. cond. water at 0oC,      [W m-1 K-1]
      'k_ice', 2.10, ...         % therm. cond. ice at 0oC,        [W m-1 K-1]
      'Rv', 461.0, ...           % gas constant water vapor        [J kg-1 K-1]
      'emiss', 0.98, ...         % surface emissivity,             [-]
      'SB', 5.6696e-8, ...       % stefan boltzman constant,       [W m-2 K-4]
      'gravity', 9.81, ...       % gravitational accel.,           [m s-2]
      'kappa', 0.4, ...          % von Karman constant,            [-]
      'one_atmos', 101300.0, ... % one atm,                        [Pa = J m-3 = kg m-1 s-2]
      'epsilon', 0.622, ...      % Rd/Rv                           [-]
      'es0', 611, ...            % Ref sat vap press in C-C Eqn    [Pa]
      'S0', 1367, ...            % Solar constant in               [W m-2]
      'N0', 0.08, ...            % Marshall-Palmer parameter       [cm-4]
      'psychro', 64.6, ...       % Psychrometric constant          [Pa K-1]
      'dalr', 9.800, ...         % dry adiabatic lapse rate        [K/km]
      'malr', 5.1, ...           % moist adiabatic lapse rate      [K/km]
      'cv_air', 1280.1, ...      % volumetric heat cap dry air     [J/m3/K]
      'cv_ice', 1927534.0, ...   % volumetric heat cap ice         [J/m3/K]
      'cv_liq', 4180000.0, ...   % volumetric heat cap liquid h20  [J/m3/K]
      'fcp', 100.0, ...          % freezing curve parameter        [K-1]
      'scale_ht', 8500.0, ...    % scale height                    [m]
      'hrsperday', 24, ...       % hours per day                   [hr]
      'roLv', 3187500.0, ...     % air density x latent heat vapor [J/m3]
      'roLs', 3613350.0, ...     % air density x latent heat subl. [J/m3]
      'roLf', 306278000.0, ...   % ice density x latent heat fusion[J/m3]
      'ro_iwe', 0.917, ...       % ice water equivalent            [dimensionless]
      'ro_wie', 1.09050, ...     % liq ice equivalent              [dimensionless]
      'fcpsq', 10000, ...        % freezing curve parameter squared[-]
      'emissSB', 5.5562080e-8 ...% emissivity of water times stefan boltzman constant [W/m2/K4]
      );

   if (nargout == 1 && nargin == 0) || (strcmp('all', varargin{1}))
      varargout{1} = constants;
   end
   for n = 1:nargin
      arg = varargin{n};
      varargout{n} = constants.(arg);
   end
end

