function [opt, args, nargs] = parseoptarg(args, validopts, defaultopt)
   %PARSEOPTARG parse optional scalar text parameter in variable argument list.
   %
   %  [OPT, ARGS, NARGS] = PARSEOPTARG(ARGS, VALIDOPTS, DEFAULTOPT)
   %
   % Description
   %  [OPT, ARGS, NARGS] = PARSEOPTARG(ARGS, VALIDOPTS, DEFAULTOPT) returns OPT,
   %  a char contained in VALIDOPTS found in ARGS, a new version of ARGS with
   %  OPT removed, and NARGS, the number of returned arguments in ARGS. If no
   %  occurences of VALIDOPTS are found in ARGS, OPT is set to DEFAULTOPT.
   %
   % PARSEOPTARG is intended to isolate a single scalar text value OPT in
   % functions using VARARGIN as the input argument, also known as a "flag".
   %
   % Example
   %
   %    function demo_function
   %    % Call the example calling_function
   %    calling_function('option2', 42, 'hello'); % 'option2' is selected
   %    calling_function(42, 'hello'); % 'option1' is selected as the default
   %    end
   %
   %    function calling_function(varargin)
   %    valid_options = {'option1', 'option2'};
   %    default_option = 'option1';
   %
   %    [selected_option, remaining_args, nargs] = parseoptarg( ...
   %       varargin, valid_options, default_option);
   %
   %    disp(['Selected option: ', selected_option]);
   %    disp(['Number of Remaining arguments: ', num2str(nargs)]);
   %    disp('Remaining arguments:');
   %    disp(remaining_args);
   %    end
   %
   % See also parseparampairs
   %
   %#codegen

   % Note: this pass the inputs in as "varargin" rather than "varargin{:}"
   if ischar(args) && isrow(args)
      args = {args};
   end

   %  PARSE INPUTS
   if nargin < 3
      defaultopt = '';
   end

   [args{1:numel(args)}] = convertStringsToChars(args{:});
   validopts = tocellstr(validopts);

   %  MAIN
   for thisarg = transpose(validopts(:))
      % Find possible char opts and remove the matching one if found.
      iopt = cellfun(@(a) ischar(a), args);
      iopt(iopt) = cellfun(@(a) strcmp(a, thisarg), args(iopt));
      opt = args(iopt);
      args = args(~iopt);
      nargs = numel(args);

      if ~isempty(opt)
         opt = opt{:}; % since optarg is a scalar text, this should work
         break
      end
   end

   %  PARSE OUTPUTS
   if isempty(opt)
      opt = defaultopt; % Initialize to default arg
   else
      if islogical(defaultopt)
         opt = true;
      end
   end
end
