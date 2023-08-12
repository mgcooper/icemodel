function [ice1, ice2] = SAVEOUTPUT(ice1, ice2, Tsfc, Qm, Qf, Qe, Qh, Qc, chi, ...
   balance, dt_sum, T, f_ice, f_liq, d_liq, d_drn, d_evp, Sc, errH, errT, iter)

% for a stripped-down run this is all that's needed:
% df_liq(:,iter)    =  d_liq;
% frac_ice(:,iter)  =  f_ice;         % fraction ice
% T_ice(:,iter)     =  T;
% T_srf(iter,1)     =  Tsfc;

% save the surface energy balance
ice1.Tsfc(iter,1)       =  Tsfc;          % surface temp
ice1.Qm(iter,1)         =  Qm;            % melt energy
ice1.Qf(iter,1)         =  Qf;            % freeze deficit
ice1.Qe(iter,1)         =  Qe;            % latent
ice1.Qh(iter,1)         =  Qh;            % sensible
ice1.Qc(iter,1)         =  Qc;            % conductive into surface
ice1.chi(iter,1)        =  chi;           % skin parameter
ice1.balance(iter,1)    =  balance;       % SEB residual
ice1.dt(iter,1)         =  dt_sum;        % check
%ice1.zD(iter,1)        =  zD;

% Save the ice column data
ice2.Tice(:,iter)       =  T;             % ice temperature
ice2.f_ice(:,iter)      =  f_ice;         % fraction ice
ice2.f_liq(:,iter)      =  f_liq;         % fraction liq
ice2.df_liq(:,iter)     =  d_liq;
ice2.df_drn(:,iter)     =  d_drn;
ice2.df_evp(:,iter)     =  d_evp;
ice2.Sc(:,iter)         =  Sc;            % source term
ice2.errH(:,iter)       =  errH;          % enthalpy error
ice2.errT(:,iter)       =  errT;          % temperature error

% save the diagnostic data
% diags.Tflag(iter,1)     =  Tflag;
% diags.LCflag(iter,1)    =  LCflag(1);

% for a stripped-down run this is all that's needed:
% T_sfc    (iter,1)    =  Tsfc;
% T_ice    (:,iter)    =  T;
% frac_ice (:,iter)    =  f_ice;
% frac_liq (:,iter)    =  f_liq;
% df_liq   (:,iter)    =  d_liq;