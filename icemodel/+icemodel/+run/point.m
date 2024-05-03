function [ice1, ice2, met, opts] = point(kwargs)

   arguments (Input)
      kwargs.saveflag (1, 1) logical = false
      kwargs.sitename (1, :) string {mustBeMember(kwargs.sitename, ...
         ["sector", "behar", "ak4", "slv1", "slv2", "upperbasin", "kanm", "kanl"])} = []
      kwargs.forcings (1, :) string {mustBeMember(kwargs.forcings, ...
         ["mar", "kanm", "kanl"])} = "mar"
      kwargs.userdata (1, :) string {mustBeMember(kwargs.userdata, ...
         ["mar", "modis", "merra", "racmo", "kanm", "kanl"])} = []
      kwargs.uservars (1, :) string {mustBeMember(kwargs.uservars, ...
         ["albedo"])} = []
      kwargs.smbmodel (1, :) string {mustBeMember(kwargs.smbmodel, ...
         ["icemodel", "skinmodel"])} = "skinmodel"
      kwargs.simyears (1, :) double = []
      kwargs.gridcell (1, :) double = []
      kwargs.testname (1, :) string = []
   end
   [saveflag, sitename, forcings, userdata, uservars, smbmodel, ...
      simyears, gridcell, testname] = dealout(kwargs);

   if isempty(userdata)
      userdata = forcings;
   end

   % varargin = struct2cell(kwargs);
   % [varargin{:}] = convertStringsToChars(varargin{:});

   %% Set the model options
   opts = icemodel.setopts(smbmodel, sitename, simyears, forcings, ...
      userdata, uservars, saveflag);

   if notempty(gridcell)
      opts.metfname = {fullfile(opts.pathinput, 'met', 'sector', ...
         ['met_' int2str(gridcell) '.mat'])};
   end

   % run the model
   switch smbmodel
      case 'icemodel'
         tic; [ice1, ice2] = icemodel(opts); toc
      case 'skinmodel'
         tic; [ice1, ice2] = skinmodel(opts); toc
   end

   % load the met data and run the post processing
   if saveflag
      [ice1, ice2, met] = icemodel.loadresults(opts);
   else
      [ice1, ice2, met] = POSTPROC(ice1, ice2, opts, simyears);
   end

   %% special case of point-scale sector met run
   if strcmp(sitename, 'sector') == true

      figure
      tiledlayout('flow')

      for var = ["Tsfc", "runoff", "melt"]
         nexttile
         plot(ice1.Time, ice1.(var))
         ylabel(var)
      end
      return
   end

   %% prep the output for plotting
   [Runoff, Discharge, Catchment] = prepRunoff(opts, ice1);
   AblationHourly = prepAblation(opts, ice1, 'hourly');
   AblationDaily = prepAblation(opts, ice1, 'daily');

   % this controls the time period over which ablation and runoff are plotted
   t1 = datetime(simyears(1),6,1,0,0,0,'TimeZone','UTC');
   t2 = datetime(simyears(1),9,1,0,0,0,'TimeZone','UTC');

   %% Plot runoff
   [h1, data] = plotRunoff(Runoff, Discharge, Catchment, ...
      'sitename', sitename, 'userdata', userdata, 'forcings', forcings, ...
      't1', t1, 't2', t2, 'refstart', false, 'smbmodel', smbmodel);

   % plot ablation
   h2 = plotPromice(AblationDaily,'refstart',t1);
   h3 = plotPromice(AblationHourly,'refstart',t1);

   % Repeat with July
   t1 = datetime(simyears(1),7,1,0,0,0,'TimeZone','UTC');
   h2 = plotPromice(AblationDaily,'refstart',t1);
   h3 = plotPromice(AblationHourly,'refstart',t1);

   %%

   dz_therm = opts.dz_thermal;
   z0_therm = opts.z0_thermal;
   [dz, delz] = CVMESH(z0_therm, dz_therm);

   %%
   % % plot the energy balance
   % plotEnbal(ice1, met);
   %
   % diffs = (data{end, 2:end} - data.ADCP(end)) ./ data.ADCP(end);
   % diffs = array2table(diffs, 'VariableNames', {'Racmo', 'Mar', 'Merra', 'IceModel'})
   % [ice1.runoff(end) ice1.melt(end)]
   %
   % figure
   % plot(ice1.balance)
   % ylabel('balance')
   %
   % figure
   % plot(ice1.chf)
   % ylabel('Q_c')
   %
   % figure
   % histogram(ice1.Tice_numiter, 'Normalization', 'probability')
   % ylabel('ICEENBAL numiter')
   %
   % fprintf('number successful Tsfc: %.0f\n', sum(ice1.Tsfc_converged))
   % fprintf('percent icenbal < 5 iterations: %.2f\n', ...
   %    100 * sum(ice1.Tice_numiter < 5) / numel(ice1.Tice_numiter))
   %
   %
   % % figure;
   % % plot(ice1.Time, ice1.melt); hold on
   % % plot(met.Time, cumsum(met.melt))
   % % plot(ice1.Time, ice1.runoff);
   % % legend('total melt', 'met melt', 'runoff')
   %
   % % zD = ice1.zD;
   % % [min(zD(zD > 0)) mean(zD(zD > 0)) max(zD(zD > 0))]
   % %
   % % figure; plot(ice1.Time, zD); hold on
   % % horzline(opts.dz_thermal)
   %
   %
   % %%
   % % icemodel.diags.stability(met.Time, met.tair, ice1.tsfc, met.wspd, ...
   % %       opts.z_0, opts.z_obs, opts.z_wind)
   %
   % %%
   %
   % % [Lv, ro_liq] = icemodel.physicalConstant('Lv', 'ro_liq');
   % % evap = ice1.lhf ./ (Lv * ro_liq) .* opts.dt ./ opts.dz_thermal;
   % %
   % % figure;
   % % plot(ice1.Time, ice1.evap); hold on;
   % % plot(ice1.Time, cumsum(evap, 'omitnan'));
   %
   % % postive lhf means energy into the surface
end
