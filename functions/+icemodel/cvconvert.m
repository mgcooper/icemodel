function varargout = cvconvert(varargin)
   %CVCONVERT Convert between control volume properties.
   % 
   % Convert constituent volume properties among several dimensions including
   % mass, volume fraction, bulk density, total density, and volume in each CV.
   %
   % Syntax:
   % [out1, ..., outN] = cvconvert(FROM, TO, CV, CONSTS, in1, ..., inN) Converts
   % numeric values in1, ..., inN from their value in FROM to their value in TO
   % using the control volumes in CV and the rquired physical constants in 
   % CONSTS.
   %
   % Inputs: 
   % convertFrom - String indicating the dimension of the input quantities
   %               ('mass', 'volumefraction', 'bulkdensity', 'totaldensity', or
   %               'volume').
   % convertTo - String indicating the dimension of the output quantities
   %             ('mass', 'volumefraction', 'bulkdensity', 'totaldensity', or
   %             'volume').
   % dz - Scalar or vector indicating the CV size(s). CV - Vector or matrix of
   % constants (i.e., densities, heat capacities) needed for conversions. 
   % varargin - One or more vectors or matrices, each corresponding to a
   % constituent, with elements representing the constituent's quantity in the
   % dimension specified by convertFrom.
   %
   % Outputs: 
   % varargout - One or more vectors or matrices, each corresponding to a
   % constituent, with elements representing the constituent's quantity in the
   % dimension specified by convertTo.
   %
   % Note: 'volumefraction' and 'bulkdensity' require the total volume or mass,
   % which is assumed to be the sum of the inputs.
   %
   % Example: 
   % [m_liq, m_ice] = icemodel.cvconvert( ...
   %    'volumefraction', 'mass', dz, [ro_liq, ro_ice], f_liq, f_ice); 
   % This converts the volume fractions of liquid
   % and ice into their corresponding masses.
   %
   % See also: icemodel, cvpropertylist

   % Parse input arguments
   assert(nargin >= 5, 'Insufficient number of input arguments.');

   convertFrom = varargin{1};
   convertTo = varargin{2};
   CV = varargin{3};
   constants = varargin{4};
   values = varargin(5:end);

   % Ensure the size of the CV matches the size of the fractions
   if ~isscalar(CV)
      assert(all(size(CV) == size(values{1})), ...
         'Size of CV must match size of input fractions.');
   end

   % Perform conversion
   switch convertFrom
      case 'volumefraction' % fk
         switch convertTo
            case 'mass' % mk = fk * rok * vt
               varargout = cell(size(values));
               for i = 1:numel(values)
                  % Volumetric fraction to mass conversion
                  varargout{i} = values{i} .* constants(i) .* CV;
               end
            case 'bulkdensity' % gk = fk * rok
               varargout = cell(size(values));
               for i = 1:numel(values)
                  % Volumetric fraction to bulk density conversion
                  varargout{i} = values{i} .* constants(i);
               end
            case 'volume' % vk = fk * vt
               varargout = cell(size(values));
               for i = 1:numel(values)
                  % Volumetric fraction to volume conversion
                  varargout{i} = values{i} .* CV;
               end
            otherwise
               error('Unsupported conversion type.');
         end
      case 'bulkdensity' % gk
         switch convertTo
            case 'totaldensity' % sum of gk
               varargout = cell(1);
               % Compute total density
               totalDensity = 0;
               for i = 1:numel(values)
                  totalDensity = totalDensity + values{i};
               end
               varargout{1} = totalDensity;
            case 'volumefraction' % fk = gk / rok
               varargout = cell(size(values));
               for i = 1:numel(values)
                  % Bulk density to volumetric fraction conversion
                  varargout{i} = values{i} ./ constants(i);
               end
            case 'volume' % vk = gk / rok * vt
               varargout = cell(size(values));
               for i = 1:numel(values)
                  % Bulk density to volume conversion
                  varargout{i} = values{i} .* CV ./ constants(i);
               end
            case 'mass' % mk = gk * vt
               varargout = cell(size(values));
               for i = 1:numel(values)
                  % Bulk density to mass conversion
                  varargout{i} = values{i} .* CV;
               end

            case 'massfraction' % fmk = gk * sum(gk) = mk / mt
               varargout = cell(size(values));
               for i = 1:numel(values)
                  % Bulk density to mass fraction conversion
                  varargout{i} = values{i} ./ sum(values{i});
               end
            otherwise
               error('Unsupported conversion type.');
         end
      case 'massfraction' % mk / (vt * sum of gk)
         switch convertTo
            case 'mass'
               error('conversion from massfraction to mass is not supported')
            otherwise
               error('Unsupported conversion type.');
         end
      case 'mass' % mk
         switch convertTo
            case 'volumefraction' % fk = mk / (rok * vt)
               varargout = cell(size(values));
               for i = 1:numel(values)
                  % Mass to volumetric fraction conversion
                  varargout{i} = values{i} ./ (constants(i) .* CV);
               end
            case 'bulkdensity' % gk = mk / vt
               varargout = cell(size(values));
               for i = 1:numel(values)
                  % Mass to bulk density conversion
                  varargout{i} = values{i} ./ CV;
               end
            case 'volume' % vk = mk / rok
               varargout = cell(size(values));
               for i = 1:numel(values)
                  % Mass to volume conversion
                  varargout{i} = values{i} ./ constants(i);
               end
            case 'massfraction' % fmk = mk / sum(mk) = mk / mt
               varargout = cell(size(values));
               for i = 1:numel(values)
                  % Mass to mass fraction conversion
                  varargout{i} = values{i} ./ sum([values{:}]);
               end
            otherwise
               error('Unsupported conversion type.');
         end
      otherwise
         error('Unsupported conversion type.');
   end
end
