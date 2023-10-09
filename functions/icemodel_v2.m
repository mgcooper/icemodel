function [ice1, ice2, met, opts] = icemodel(opts) %#codegen
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

%% model initialization

% Load the physical constants
[  cv_air, cv_liq, cv_ice, emiss, SB, epsilon, fcp, k_liq, Lf, Ls, Lv, ro_air, ...
   ro_ice, ro_liq, roLf, roLs, roLv, Rv, Tf] = icemodel.physicalConstant( ...
   'cv_air','cv_liq','cv_ice','emiss', 'SB','epsilon','fcp','k_liq','Lf', ...
   'Ls','Lv', 'ro_air','ro_ice','ro_liq','roLf','roLs','roLv', 'Rv','Tf');
TINY = 1e-8;

% Load the forcing data
[met, opts] = METINIT(opts, 1);

% Initialize the ice column
[f_ice, f_liq, T, TL, TH, flmin, flmax, cp_sno, ~, dz, fn, delz, grid_therm, ...
   dz_therm, dz_spect, JJ_therm, JJ_spect, ~, Sp, wcoef, scoef, ro_sno, ...
   ro_iwe, ro_wie, xTsfc, xf_liq, roL, Qc, f_min, fopts, liqflag, ...
   ice1, ice2] = ICEINIT(opts, met);

% Initialize the extinction coefficients
[total_solar, grid_spect, z_walls, spect_lower, spect_upper, ...
   solardwavl] = EXTCOEFSINIT(opts, dz_spect, JJ_spect, ro_ice);

% Initialize timestepping
[iter, subiter, itime, maxiter, maxsubiter, dt, dt_min, ...
   dt_max, dt_new] = INITTIMESTEPS(opts, met);

%% Start iterations over years

for MM = 1:opts.annual_loops

   while iter <= maxiter

      % Initialize substep
      [dt_sum, subfail, dt_flag, OK, d_liq, d_drn, d_evp] = INITSUBSTEP(f_liq);
      
      % Load the met data for this timestep
      [Tair, rh, wspd, Qsi, Qli, Pa, albedo, De, ea] = LOADMETDATA( ...
         met, iter, wcoef, liqflag);

      % update the upper layer diffusion length scale
      % if sqrt(k_eff(1)*dt/(ro_sno(1)*cp_sno(1))) > dz(1)
      %    pause;
      % end

      % Update the subsurface solar radiation source term
      [Sc, chi] = UPDATEEXTCOEFS(grid_spect, JJ_spect, dz_spect, ...
         grid_therm, JJ_therm, dz_therm, z_walls, dz, ro_sno, Qsi, ...
         total_solar, spect_lower, spect_upper, albedo, solardwavl);

      % could use this to get a better estimate of Tsfc if convergence fails
      % xTsfc = T(1) + SEB/(k_eff(1)/(dz/2)); % T1 + SEB/a1

      while OK == false || dt_sum < dt

         % Surface temperature
         Tsfc = SEBSOLVE(Tair, Qsi, Qli, ea, albedo, De, Pa, wspd, cv_air, ...
            emiss, SB, Tf, Qc, xTsfc, chi, roL, scoef, fopts, liqflag);

         % Subsurface energy balance
         [T, errH, errT, f_ice, f_liq, OK] = ICEENBAL(T, f_ice, f_liq, ...
            k_liq, cv_ice, cv_liq, ro_ice, ro_liq, ro_sno, cp_sno, Ls, ...
            Lf, roLf, Rv, Tf, dz, delz, fn, dt_new, JJ_therm, Tsfc, Sc, ...
            fcp, TL, TH, flmin, flmax, ro_iwe, ro_wie);

         % Error message (slows down the code a lot)
         % fprintf('iter = %d (%.2f%%), dt = %.0f, success = %s\n',          ...
         %    iter,100*iter/maxiter,dt_new,mat2str(OK))

         % Phase boundary overshoot, decrease the timestep and start over
         if OK == false && subfail < maxsubiter
            [subfail, subiter, dt_new, T, Tsfc, f_ice, f_liq] = RESETSUBSTEP( ...
               xT, xTsfc, xf_ice, xf_liq, dt_max, subiter, maxsubiter, subfail);
            continue
         else

            % Update effective thermal conductivity
            k_eff = GETGAMMA(T,f_liq,f_ice,ro_ice,k_liq,Ls,Rv,Tf);

            % Compute surface fluxes
            [Qe, Qh, Qc, Qm, Qf, ea, balance] = SEBFLUX(T, Tair, Tsfc, ...
               Qsi, Qli, albedo, wspd, rh, Pa, Tf, k_eff, cv_air, roL, ...
               emiss, SB, epsilon, scoef, De, dz, liqflag, chi);

            % Compute melt freeze
            [T, f_ice, f_liq, d_liq, d_evp, d_drn] = ICEMF(T, f_ice, ...
               f_liq, ro_ice, ro_liq, cv_ice, cv_liq, Lf, Ls, Lv, ...
               Tf, TL, fcp, xf_liq, Sc, Sp, JJ_therm, f_min, fopts, ...
               dz_therm, dt_new, Qe, liqflag, ro_iwe, d_liq, d_drn, ...
               d_evp, flmin, iter);

            % Update density, heat capacity, and substep time 
            [ro_sno, cp_sno, liqflag, roL, xT, xTsfc, xf_liq, xf_ice, ...
               dt_sum, itime, dt_new, dt_flag] = UPDATESUBSTEP(f_ice, ...
               f_liq, ro_ice, ro_liq, ro_air, cv_ice, cv_liq, T, Tsfc, ...
               dt, dt_sum, itime, dt_new, roLv, roLs, dt_min, TINY);
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
   [iter, subiter, itime] = INITTIMESTEPS(opts, met);
end
% Post process
[ice1,ice2,met] = POSTPROC(ice1,ice2,met,opts);
