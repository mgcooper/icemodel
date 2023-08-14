function [ice1, ice2, opts] = icemodel(opts) %#codegen
% ICEMODEL Simulate the phase change process in melting glacier ice.
%
% This function models the phase change process in melting glacier ice. It uses
% iterative processes to update the temperature, liquid and ice fraction at
% each time step. This model considers both the surface and subsurface energy
% balance to simulate the phase change.
%
% Syntax:
% [ice1, ice2, met, opts] = ICEMODEL(opts)
%
% Inputs:
% opts - A structure containing model options and parameters. Defined by the
%        icemodel_opts function. The fields include:
%        * sitename         - Site name for the model simulation
%        * simmodel         - Simulation model identifier
%        * simyears         - Years for which the simulation is done
%        * forcings         - Type of forcing data used
%        * userdata         - User-defined data type
%        * uservars         - User-defined variables
%        * savedata         - Flag indicating if data should be saved
%        * testname         - Name of the test (default: 'none')
%        * testpoint        - Test point if applicable (default: 'none')
%        * ... (other parameters related to the model configuration)
%
% Outputs:
% ice1  - 1-dimensional data storing variables defined at the ice surface or
%         near-surface atmosphere. Contains one value per timestep.
% ice2  - 2-dimensional data storing variables defined on the subsurface ice
%         column control volume mesh. Contains one column per timestep.
% met   - Structure containing the meteorological data used for the simulation.
% opts  - Updated structure containing model options and parameters.
%
% See also: ICEMODEL_OPTS

%% INITIALIZE THE MODEL

% Load the physical constants
[  cv_air, cv_liq, cv_ice, emiss, SB, epsilon, fcp, k_liq, Lf, ...
   Ls, Lv, ro_air, ro_ice, ro_liq, roLf, roLs, roLv, Rv, Tf] = ...
   icemodel.physicalConstant( ...
   'cv_air','cv_liq','cv_ice','emiss', 'SB','epsilon','fcp','k_liq','Lf', ...
   'Ls','Lv', 'ro_air','ro_ice','ro_liq','roLf','roLs','roLv', 'Rv','Tf');
TINY = 1e-8;

% Load the forcing data
[opts, tair, swd, lwd, albedo, wspd, rh, psfc, De, Time] = METINIT(opts, 1);

% Initialize the ice column
[f_ice, f_liq, T, TL, TH, flmin, flmax, cp_sno, ~, dz, fn, delz, grid_therm, ...
   dz_therm, dz_spect, JJ_therm, JJ_spect, ~, Sp, scoef, ro_sno, ro_iwe, ...
   ro_wie, xTsfc, xf_liq, roL, Qc, f_min, fopts, liqflag, ice1, ice2] = ...
   ICEINIT(opts, tair(1));

% Initialize the extinction coefficients
[total_solar, grid_spect, z_walls, spect_lower, spect_upper, solardwavl] = ...
   EXTCOEFSINIT(opts, dz_spect, JJ_spect, ro_ice);

% Initialize timestepping
[iter, subiter, itime, maxiter, maxsubiter, dt, dt_min, dt_max, dt_new] = ...
   INITTIMESTEPS(opts, Time);

%% START ITERATIONS OVER YEARS

for MM = 1:opts.annual_loops

   while iter <= maxiter

      % INITIALIZE NEW SUBSTEP
      [dt_sum, subfail, dt_flag, OK, d_liq, d_drn, d_evp] = ...
         INITSUBSTEP(f_liq);

      % update the upper layer diffusion length scale
      % if sqrt(k_eff(1)*dt/(ro_sno(1)*cp_sno(1))) > dz(1)
      %    pause;
      % end

      % SUBSURFACE SOLAR RADIATION SOURCE-TERM
      if swd(iter) > 0
         [Sc, chi] = UPDATEEXTCOEFS(swd(iter), albedo(iter), grid_spect, ...
            JJ_spect, dz_spect, grid_therm, JJ_therm, dz_therm, dz, z_walls, ...
            ro_sno, total_solar, spect_lower, spect_upper, solardwavl);
      else
         Sc = 0.0.*dz;
         chi = 0.0;
      end

      % [Qsi*(1-albedo) sum(Sc.*dz)]                        % total
      % [chi*Qsi*(1-albedo) Sc(1).*dz(1)]                   % top node
      % [(1-chi)*Qsi*(1-albedo) sum(Sc(2:end).*dz(2:end))] 	% interior nodes

      % if the Sc(1)=0 method is used, then this should hold:
      % Qsi*(1-albedo) - ( chi*Qsi*(1-albedo) + sum(Sc.*dz) ) = 0
      % (1-chi)*Qsi*(1.0-albedo) - Sc(1)*dz(1) = 0

      % if using the Qseb method, the total absorbed radiation should equal the
      % portion allocated to the surface plus the sum of the subsurface absorption
      % [Qsi*(1-albedo) sum(Sc.*dz)+Qseb]                   % total
      % [chi*Qsi*(1.0-albedo) Qseb]                         % SEB
      % [(1-chi)*Qsi*(1.0-albedo) sum(Sc.*dz)]              % interior

      % so we can either pass Qseb or chi out of extcoefs

      % SURFACE TEMPERATURE
      ea = VAPPRESS(rh(iter), tair(iter), liqflag);
      Tsfc = fsearchzero(@(Tsfc) ...
         chi * (1.0-albedo(iter)) * swd(iter) + emiss * lwd(iter) + Qc - ...
         emiss * SB * Tsfc^4 + cv_air * De(iter) * ...
         STABLEFN(tair(iter), Tsfc, wspd(iter), scoef) * ...
         (tair(iter)-Tsfc) + roL * De(iter) * 0.622 / psfc(iter) * ...
         STABLEFN(tair(iter), Tsfc, wspd(iter), scoef) * ...
         (ea - VAPOR(Tsfc, Tf, liqflag)), ...
         xTsfc, xTsfc-50, xTsfc+50, tair(iter), fopts.TolX);

      while OK == false || dt_sum < dt

         % SUBSURFACE ENERGY BALANCE
         [T, errH, errT, f_ice, f_liq, OK] = ICEENBAL(T, f_ice, f_liq, ...
            k_liq, cv_ice, cv_liq, ro_ice, ro_liq, ro_sno, cp_sno, Ls, ...
            Lf, roLf, Rv, Tf, dz, delz, fn, dt_new, JJ_therm, Tsfc, Sc, ...
            fcp, TL, TH, flmin, flmax, ro_iwe, ro_wie);

         % ERROR MESSAGE (SLOWS DOWN THE CODE A LOT)
         % fprintf('iter = %d (%.2f%%), dt = %.0f, success = %s\n',          ...
         %    iter,100*iter/maxiter,dt_new,mat2str(OK))

         % PHASE BOUNDARY OVERSHOOT, DECREASE THE TIME STEP AND START OVER
         if not(OK) && subfail < maxsubiter
            [subfail, subiter, dt_new, T, Tsfc, f_ice, f_liq] = ...
               RESETSUBSTEP( xT, xTsfc, xf_ice, xf_liq, dt_max, subiter, ...
               maxsubiter, subfail);
            continue
         else

            % UPDATE SURFACE FLUXES
            k_eff = GETGAMMA(T, f_liq, f_ice, ro_ice, k_liq, Ls, Rv, Tf);
            [Qe, Qh, Qc, Qm, Qf, balance] = SEBFLUX(T, Tsfc, tair(iter), ...
               swd(iter), lwd(iter), albedo(iter), wspd(iter), psfc(iter), ...
               De(iter), ea, Tf, k_eff, cv_air, roL, emiss, SB, epsilon, ...
               scoef, dz, liqflag, chi);

            % COMPUTE MELT FREEZE
            [T, f_ice, f_liq, d_liq, d_evp, d_drn] = ICEMF(T, f_ice, ...
               f_liq, ro_ice, ro_liq, cv_ice, cv_liq, Lf, Ls, Lv, ...
               Tf, TL, fcp, xf_liq, Sc, Sp, JJ_therm, f_min, fopts, ...
               dz_therm, dt_new, Qe, liqflag, ro_iwe, d_liq, d_drn, ...
               d_evp, flmin);

            % UPDATE DENSITY, HEAT CAPACITY, AND SUBSTEP TIME
            [ro_sno, cp_sno, liqflag, roL, xT, xTsfc, xf_liq, xf_ice, ...
               dt_sum, itime, dt_new, dt_flag] = ...
               UPDATESUBSTEP(f_ice, f_liq, ro_ice, ro_liq, ro_air, cv_ice, ...
               cv_liq, T, Tsfc, dt, dt_sum, itime, dt_new, roLv, roLs, ...
               dt_min, TINY);
            
            % zD = sqrt(k_eff(1)*dt/(ro_sno(1)*cp_sno(1)));
            % if zD > dz(1)
            %
            % end
         end
      end
      
      % Save output
      if MM==opts.annual_loops
         [ice1, ice2] = SAVEOUTPUT(ice1, ice2, Tsfc, Qm, Qf, Qe, Qh, Qc, ...
            chi, balance, dt_sum, T, f_ice, f_liq, d_liq, d_drn, d_evp, ...
            Sc, errH, errT, iter);
      end
      
      % Move to the next timestep
      [iter, subiter, dt_new] = NEXTSTEP(iter, subiter, dt_flag, dt_max, OK);
   end
   
   % Restart the model (spinup)
   [iter, subiter, itime] = INITTIMESTEPS(opts, Time);
end

% Post process
[ice1, ice2] = POSTPROC(ice1, ice2, opts, swd, lwd, albedo, Time);
