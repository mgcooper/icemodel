function sample = sampleMachineState(kwargs)
   %SAMPLEMACHINESTATE Sample the machine conditions that disturb a timing.
   %
   %  sample = icemodel.test.helpers.sampleMachineState()
   %  sample = icemodel.test.helpers.sampleMachineState(own_pid=1234)
   %
   % run_perf_suite and build_perf_baseline call this at run start, after
   % each case, and at run end. summarizeMachineState reduces the samples
   % to the attestation that the artifact records.
   %
   % Fields of SAMPLE
   %  foreign_matlab_processes  Count of MATLAB processes that are not this
   %                            process and do not descend from it. NaN
   %                            when the process list cannot be read.
   %  load_average_1min         One-minute load average on macOS and Linux.
   %                            Windows has no load average, so the field
   %                            holds the processor queue length, the count
   %                            of threads waiting for a processor, which
   %                            measures the same contention. NaN when it
   %                            cannot be read.
   %  ac_power                  True when the machine draws AC power. macOS
   %                            reads pmset; Windows reads Win32_Battery,
   %                            and a machine without a battery reports
   %                            true; Linux reads the AC online flag under
   %                            /sys/class/power_supply, and a machine
   %                            without a battery reports true.
   %  sampled_utc               Sample time.
   %  probe_errors              Text of every probe that failed, so a NaN
   %                            field can be explained after the fact.
   %
   % Name-value
   %  command_runner  Runs each probe command; returns [status, text] like
   %                  system. Tests inject probe text through it.
   %  own_pid         Process id of this MATLAB. Defaults to feature getpid.
   %  platform        "mac", "linux", or "windows". Defaults to the running
   %                  platform. Tests exercise each branch through it.
   %  power_supply_dir  Root of the Linux power-supply tree.
   %
   % See also: icemodel.test.helpers.summarizeMachineState,
   %  icemodel.test.helpers.perfMeasurementQuality

   arguments
      kwargs.command_runner (1, 1) function_handle = @system
      kwargs.own_pid (1, 1) double = feature('getpid')
      kwargs.platform (1, 1) string ...
         {mustBeMember(kwargs.platform, ["mac", "linux", "windows"])} ...
         = runningPlatform()
      kwargs.power_supply_dir (1, 1) string = "/sys/class/power_supply"
   end

   probe_errors = strings(0, 1);
   platform = kwargs.platform;

   % A MATLAB process is one whose command names the MATLAB binary under a
   % platform bin folder. The MCP servers and the matlab launcher script
   % mention MATLAB in their paths but are not MATLAB, so the pattern
   % requires the binary itself. Both probes print one "pid ppid command"
   % line per process, so one parser serves every platform.
   if platform == "windows"
      process_command = "powershell -NoProfile -Command """ + ...
         "Get-CimInstance Win32_Process | ForEach-Object { " + ...
         "'{0} {1} {2}' -f $_.ProcessId, $_.ParentProcessId, $_.CommandLine }""";
   else
      process_command = "ps -axo pid,ppid,command";
   end
   [status, text] = kwargs.command_runner(process_command);
   foreign = NaN;
   if status == 0
      foreign = countForeignMatlab(string(text), kwargs.own_pid);
   else
      probe_errors(end + 1, 1) = "process list: " + string(strtrim(text));
   end

   % The load probe differs by platform, so resolve it in one local function.
   [load_average, load_error] = readLoad(kwargs.command_runner, platform);
   if strlength(load_error) > 0
      probe_errors(end + 1, 1) = load_error;
   end

   % The power probe differs by platform, so resolve it in one local function.
   [ac_power, power_error] = readAcPower(kwargs.command_runner, ...
      platform, kwargs.power_supply_dir);
   if strlength(power_error) > 0
      probe_errors(end + 1, 1) = power_error;
   end

   sample = struct( ...
      'foreign_matlab_processes', foreign, ...
      'load_average_1min', load_average, ...
      'ac_power', ac_power, ...
      'sampled_utc', datetime('now', 'TimeZone', 'UTC'), ...
      'probe_errors', probe_errors);
end

function n_foreign = countForeignMatlab(ps_text, own_pid)
   %COUNTFOREIGNMATLAB Count MATLAB processes outside this process tree.
   lines = splitlines(ps_text);
   lines = lines(strlength(strtrim(lines)) > 0);
   % Each line is "pid ppid command". The header line has no numeric pid
   % and drops out of the parse. A system process can print no command at
   % all; its record must still enter the parent chain, so the command
   % group may be empty.
   tokens = regexp(lines, '^\s*(\d+)\s+(\d+)\s*(.*)$', 'tokens', 'once');
   keep = ~cellfun(@isempty, tokens);
   tokens = vertcat(tokens{keep});
   if isempty(tokens)
      n_foreign = 0;
      return
   end
   pids = str2double(tokens(:, 1));
   ppids = str2double(tokens(:, 2));
   commands = string(tokens(:, 3));
   % Windows paths are case-insensitive, and Win32_Process can report the
   % same binary as matlab.exe or MATLAB.exe, so the match ignores case on
   % every platform; no other binary shares this path shape.
   is_matlab = ~cellfun(@isempty, regexpi(commands, ...
      '[/\\]bin[/\\](maca64|maci64|glnxa64|win64)[/\\]MATLAB(\.exe)?"?(\s|$)', ...
      'once'));

   % Walk each MATLAB process up its parent chain. One that reaches this
   % process is a subprocess this run started; every other one is foreign.
   n_foreign = 0;
   for k = find(is_matlab).'
      pid = pids(k);
      own = false;
      visited = 0;
      while pid > 0 && visited < numel(pids)
         if pid == own_pid
            own = true;
            break
         end
         parent = ppids(pids == pid);
         if isempty(parent)
            break
         end
         pid = parent(1);
         visited = visited + 1;
      end
      n_foreign = n_foreign + ~own;
   end
end

function platform = runningPlatform()
   %RUNNINGPLATFORM Name the platform this MATLAB runs on.
   if ismac()
      platform = "mac";
   elseif ispc()
      platform = "windows";
   else
      platform = "linux";
   end
end

function [load_value, probe_error] = readLoad(command_runner, platform)
   %READLOAD Read the one-minute load average, or its Windows stand-in.
   probe_error = "";
   load_value = NaN;
   if platform == "windows"
      % ProcessorQueueLength counts threads waiting for a processor, the
      % closest Windows counter to a load average.
      [status, text] = command_runner("powershell -NoProfile -Command """ + ...
         "(Get-CimInstance Win32_PerfFormattedData_PerfOS_System)" + ...
         ".ProcessorQueueLength""");
      value = str2double(strtrim(string(text)));
      if status == 0 && isfinite(value)
         load_value = value;
      else
         probe_error = "processor queue: " + string(strtrim(text));
      end
      return
   end
   % uptime prints "load averages: a b c" on macOS and "load average: a, b,
   % c" on Linux; the first number is the one-minute value on both.
   [status, text] = command_runner("uptime");
   token = regexp(string(text), 'load averages?:\s*([0-9.]+)', 'tokens', 'once');
   if status == 0 && ~isempty(token)
      load_value = str2double(token{1});
   else
      probe_error = "uptime: " + string(strtrim(text));
   end
end

function [ac_power, probe_error] = readAcPower(command_runner, platform, ...
      power_supply_dir)
   %READACPOWER Read whether the machine draws AC power.
   probe_error = "";
   if platform == "mac"
      [status, text] = command_runner("pmset -g batt");
      if status ~= 0
         ac_power = false;
         probe_error = "pmset: " + string(strtrim(text));
         return
      end
      ac_power = contains(string(text), "AC Power");
      return
   end
   if platform == "windows"
      % Win32_Battery is absent on a desktop, which draws mains power.
      % BatteryStatus 2, 3, 6, 7, 8, and 9 are the states with AC connected.
      [status, text] = command_runner("powershell -NoProfile -Command """ + ...
         "$b = Get-CimInstance Win32_Battery; " + ...
         "if ($b) { $b.BatteryStatus } else { 'none' }""");
      values = strtrim(splitlines(strtrim(string(text))));
      values = values(strlength(values) > 0);
      if status ~= 0 || isempty(values)
         ac_power = false;
         probe_error = "battery: " + string(strtrim(text));
         return
      end
      % A machine with several batteries prints one status per line, and
      % every battery must report an AC-connected state.
      ac_power = all(values == "none") ...
         || all(ismember(str2double(values), [2 3 6 7 8 9]));
      return
   end

   % A machine without a battery entry draws mains power by construction.
   batteries = dir(fullfile(power_supply_dir, "BAT*"));
   if isempty(batteries)
      ac_power = true;
      return
   end
   adapters = dir(fullfile(power_supply_dir, "AC*"));
   ac_power = false;
   for k = 1:numel(adapters)
      online_file = fullfile(adapters(k).folder, adapters(k).name, "online");
      if isfile(online_file)
         ac_power = ac_power || strtrim(string(fileread(online_file))) == "1";
      end
   end
end
