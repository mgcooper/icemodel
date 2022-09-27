%--------------------------------------------------------------------------
%   UPDATE STATE VARIABLES AFTER MELT/FREEZE, LAYER COMBINATION, AND
%   GETNEWT ARE FINISHED
%--------------------------------------------------------------------------

% temporarily keeping this as _old. thenew one came from H1d_UPDATESTATE

% compute new ice, liquid water, and air volumetric fractions, new heat
% conductivity, specific heat capacity, and partial densities

function [  f_ice,                                                      ...
            f_liq,                                                      ...
            f_air,                                                      ...
            Cp_sno,                                                     ...
            xk_eff,                                                     ...
            xk_vap,                                                     ...
            ro_sno,                                                     ...
            gam_ice,                                                    ...
            aP0_old]    =   UPDATESTATE(f_ice,f_liq,f_air,T_old,Tfp,    ...
                            Cp_liq,Cp_ice,Cp_sno,xLs,xk_liq,Rv,ro_sno,  ...
                            ro_ice,ro_liq,h_tot,dt,opts)
%-------------------------------------------------------------------------- 

% don't delete commented-out statements

if opts.skin_model == true

% Update the ice matrix thermal conductivity [W m-1 K-1]
% frac_air/liq/ice/Cp_snow/ro_snow remain constant at their initial values
    xk_sno          =   GETKTHERMAL(T_old,f_ice,ro_ice,ro_sno,opts);
    xk_vap          =   0.0.*f_ice;
    xk_eff          =   xk_sno;

else
    
% update the volumetric fractions (theta in Jordan 1991, eq. 1)
    f_air           =   1-f_liq-f_ice;
    
% update bulk densities (gamma in Jordan) (mass / total volume)
    gam_ice         =   f_ice .* ro_ice;
    gam_liq         =   f_liq .* ro_liq;
  % gam_air         =   frac_air .* ro_air;
  % gam_vap         =   frac_vap .* ro_vap_sat .* frh; % for reference

% update mass (note: h_ice = volume_ice)
    mass_ice        =   f_ice .* h_tot .* ro_ice;               % [kg m-2]
    mass_liq        =   f_liq .* h_tot .* ro_liq;               % [kg m-2]
  % mass_air        =   h_air .* ro_air;                        % [kg m-2]

% mass liquid fraction (air mass is negligible) (Eq. 
%   f_liq           =   mass_liq ./ (mass_liq + mass_ice);      % [-]
%   f_liq           =   gam_liq  ./ (gam_liq  + gam_ice);       % also true

% ice matrix ('snow') density (wet) (gamma_w = gamma_i + gamma_l in Jordan)
    ro_sno         =   (mass_ice + mass_liq) ./ h_tot;

% In Jordan, bt = bl+bi+bd, but bw = bl+bi where b=gamma. So, if there are
% no dry solids, then bt=bw, meaning ro_snow = gamma_w = gamma_t
% i.e., her bt (nodal total bulk density) is my ro_snow
    
% combined specific heat capacity (ct in Jordan, Eq. 55)    % [J m-3 K-1]
    Cp_sno          =   (gam_ice.*Cp_ice+gam_liq.*Cp_liq)./ro_sno;
%   Cp_sno          =   (gam_ice.*Cp_ice+gam_liq.*Cp_liq)./gam_wat;    
%   Cw_sno          =   f_liq.*Cp_liq + (1-f_liq).*Cp_ice;  % also true

% Jordan
    d_ro_vap_dT =   VAPOR_HEAT(T_old,Tfp,Rv,xLs,f_liq,f_air);
    aP0_old     =   (ro_sno.*Cp_sno+xLs.*f_air.*d_ro_vap_dT).*h_tot./dt;
    
% new (did not make a note, but I think this is not implemented
%     aP0_old         =   ro_sno .* Cp_sno .* dy_p ./ dt;

% thermal conductivity
    xk_sno          =   GETKTHERMAL(T_old,f_ice,ro_ice,ro_sno,opts);
    
% update vapor diffusivity based on the new values of T_old
    xk_vap          =   GETKVAPOR(T_old,xLs,Rv,Tfp);

% combine thermal conductivity (ki) and vapor diffusivity (kv)
    xk_sno          =   (f_liq .* xk_liq) + (f_ice .* xk_sno); 
    xk_eff          =   (1.0-f_air) .* xk_sno + f_air .* xk_vap;

% note: gam_ice is dry snow density, needed for kthermal expressions below
   %gam_ice(gam_ice>917) = 917;
    
end
    
% Jordan notation, where k = i for ice, l for liq, v for vap, etc.: 
% rho_k     = mass_k    / volume_k  [kg/m3] (intrinsic dens.) (e.g. ro_ice)
% theta_k   = volume_k  / volume_T  [m3/m3] (fractional vol.) (e.g. frac_ice)
% gamma_k   = mass_k    / volume_T  [kg/m3] (bulk density)  (e.g. ro_snow_dry)
% gamma_k   = theta_k*rho_k
% 1         = theta_ice + theta_air + theta_liq
% rho_t     = gamma_ice + gamma_air + gamma_liq
% gamma_w   = gamma_ice + gamma_liq (b/c gamma_air negligible) (=ro_snow)
% c_t       = specific heat (J/Kg/K) (=c_w if no soil is present)
% h_t       = specific enthalpy (J/Kg)
% h_v       = h_t * rho_t (enthalpy, not defined explicitly, [J/m3]
% c_v       = c_t * rho_t (volumetric heat cap., Cp in Clark 2021) [J/m3/K]
% c_v       = Cp_snow * ro_snow = Cp = dHt/dT in Clark, undefined in Jordan

% In Jordan: Cp_sno=ct, Cw_sno=cw, f_liq=flo
% ct = (gamma_w*(f_liq.*Cp_water + (1-f_liq).*Cp_ice)./gamma_t;   
% cw = f_liq .* Cp_water + (1 - f_liq) .* Cp_ice; 

% but for snow, they simplify b/c no dry solids term, so total density is
% snow density, etc.: (using p for rho/density for compact notation):
% pt = pw = ro_sno = gam_wat    (jordan uses pt and gam_wat)
% ct = cw = Cp_sno = Cp_wat     (jordan uses ct and cw)

% NOTE: Cp_snow is called Cs in my technical document, it is the combined
% specific heat capacity of the ice + liquid mixture
% Cp_snow*ro_snow = Cs*rho_s = Cv is the volumetric specific heat capacity,
% which isn't defined in jordan, but I use it my technical document


% this just shows them lined up for pattern recognition:

% pt=pw=ro_sno = (frac_ice.*ro_ice         + frac_liq.*ro_liq); (= gam_wat)
% ct=cw=Cp_sno = (frac_ice.*ro_ice.*Cp_ice + frac_liq.*ro_liq.*Cp_liq)./ro_sno;

% pt=pw=ro_sno = (gam_ice         + gam_liq        ); (= gam_wat)
% ct=cw=Cp_sno = (gam_ice.*Cp_ice + gam_liq.*Cp_liq)./gam_wat;
% ct=cw=Cp_sno = (1-fliq).*Cp_ice +    fliq.*Cp_liq;

% where fliq = gam_liq/gam_wat = gam_liq/ro_sno;






% % my previous definition, identical to above (Jordan in text b/w Eq. 59-60)
%     Cp_sno_1    =   f_liq .* Cp_liq + (1 - f_liq) .* Cp_ice;
% 
% % technically:
%     Cp_sno_2    =   (Cp_liq.*ro_liq.*frac_liq + ...
%                             Cp_ice.*ro_ice.*frac_ice)./ro_sno;
% 
% % and to get the total apparent heat:
%     Cp_sno_3    =   (Cp_liq.*ro_liq.*frac_liq + Cp_ice.*ro_ice.*frac_ice ...
%                         + xLs.*CkT)./ro_sno;
% 
%     Cw_sno      =   f_liq .* Cp_liq + (1 - f_liq) .* Cp_ice; 
% % In Jordan: Cp_sno=ct, Cw_sno=cw, f_liq=flo
%     ct          =   (gam_liq*(f_liq.*Cp_liq + (1-f_liq).*Cp_ice))./gam_tot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% new, to check Tubini's enthalpy condition:
% Mixture enthalpy:
%     Hliq    = Cp_water.*ro_water.*frac_liq.*(T_old-Tfp);
%     Hice    = Cp_ice.*ro_ice.*frac_ice.*(T_old-Tfp);
%     Hair    = 1004.*1.275.*frac_air.*(T_old-Tfp);
%     del_ice = updateFracIce(T_old,Tfp);
%     H       = Hliq+Hice+Hair-xLf.*ro_ice.*del_ice;      % [J m-3]
%     
%     % Jordan's total water fraction (ice + liquid)
%     frac_w  = frac_liq + ro_ice./ro_water.*frac_ice;