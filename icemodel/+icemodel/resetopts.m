function opts = resetopts(opts, varargin)
   %RESETOPTS Override existing OPTS fields by name.
   %
   %  opts = icemodel.resetopts(opts, 'dt', 900, 'solver', 3)
   %
   % Derived run settings that depend on the overridden fields are cleared
   % unless they are explicitly set in the same call. The next call to
   % icemodel.configureRun then rebuilds them from the updated OPTS.
   %
   %#codegen

   narginchk(1, inf);

   if nargin == 1
      return
   end

   if rem(nargin-1, 2) ~= 0
      error('resetopts expects name/value pairs')
   end

   names = cell(1, (nargin - 1) / 2);
   k = 0;
   for n = 1:2:nargin-1
      name = varargin{n};
      if isstring(name)
         name = char(name);
      end
      if ~isfield(opts, name)
         error('unrecognized opts field: %s', name)
      end
      k = k + 1;
      names{k} = name;
      opts.(name) = varargin{n+1};
   end

   % Keep saveopts coupled to saveflag by default unless explicitly set.
   if ismember('saveflag', names) && ~ismember('saveopts', names)
      opts.saveopts = opts.saveflag;
   end

   % Update the time-lag variable if dt changed
   if ismember('dt', names) && ~ismember('tlag', names)
      opts.tlag = 6 * 3600 / opts.dt;
   end

   % Re-apply solver-dependent coupling defaults unless explicitly overridden.
   if ismember('solver', names) && ~ismember('cpl_maxiter', names)
      if opts.solver == 2
         opts.cpl_maxiter = 1;
      else
         opts.cpl_maxiter = 100;
      end
   end

   % Clear derived run settings when their source fields change. These fields
   % will be re-populated by icemodel.configureRun.
   if ismember('pathdata', names)
      if ~ismember('pathinput', names)
         opts.pathinput = [];
      end
      if ~ismember('patheval', names)
         opts.patheval = [];
      end
      if ~ismember('pathoutput', names)
         opts.pathoutput = [];
      end
      if ~ismember('pathuserdata', names)
         opts.pathuserdata = [];
      end
   end

   if any(ismember(names, {'sitename', 'smbmodel', 'testname', 'userdata'})) ...
         && ~ismember('pathoutput', names)
      opts.pathoutput = [];
   end

   if ismember('pathinput', names) && ~ismember('pathuserdata', names)
      opts.pathuserdata = [];
   end

   if any(ismember(names, {'forcings', 'userdata', 'uservars'})) ...
         && ~ismember('casename', names)
      opts.casename = [];
   end

   if any(ismember(names, {'pathinput', 'sitename', 'forcings', ...
         'simyears', 'dt'})) && ~ismember('metfname', names)
      opts.metfname = {};
   end

   if any(ismember(names, {'sitename'})) && ~ismember('output_profile', names)
      opts.output_profile = [];
   end

   if any(ismember(names, {'sitename', 'smbmodel', 'output_profile'})) ...
         && ~ismember('vars1', names)
      opts.vars1 = {};
   end
   if any(ismember(names, {'sitename', 'smbmodel', 'output_profile'})) ...
         && ~ismember('vars2', names)
      opts.vars2 = {};
   end
end
