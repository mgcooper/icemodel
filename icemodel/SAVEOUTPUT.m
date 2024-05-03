function [ice1, ice2] = SAVEOUTPUT(iter, ice1, ice2, vars1, vars2, data1, data2)
   %SAVEOUTPUT

   % Note: the order of vars1 and vars2 (set in icemodel.setopts) must match the
   % order of data1, data2 (passed in from icemodel main). It's fine if more
   % variables exist in vars1, vars2 than there are data in data1, data2, but
   % the vars1, vars2 lists must be in the same order as the data in data1,
   % data2. Thus vars1, vars2 can be set to some base set of variables in
   % setopts, and the actual variables can be adjusted in icemodel main.
   for n = 1:numel(data1)
      ice1.(vars1{n})(iter, 1) = data1{n};
   end

   for n = 1:numel(data2)
      ice2.(vars2{n})(:, iter) = data2{n};
   end
end

% function [ice1, ice2, diags] = SAVEOUTPUT(iter, ice1, ice2, vars1, vars2, ...
%       data1, data2, diags1, diags2, diagvars1, diagvars2, diagdata1, diagdata2)
%
%    if nargin > 7
%       for n = 1:numel(diagdata1)
%          diags1.(diagvars1{n})(iter, 1) = diagdata1{n};
%       end
%       for n = 1:numel(diagdata2)
%          diags2.(diagvars2{n})(:, iter) = diagdata2{n};
%       end
%    end
% end

% % or:
% function [ice1, ice2] = SAVEOUTPUT(iter, ice1, ice2, vars1, vars2, varargin)
%
% data1 = varargin(1:numel(vars1));
% data2 = varargin(numel(vars1)+1:numel(varargin));
%
% for n = 1:numel(vars1)
%    ice1.(vars1{n})(iter, 1) = data1{n};
% end
% for n = 1:numel(vars2)
%    ice2.(vars2{n})(:, iter) = data2{n};
% end
%
%
% function [ice1, ice2] = SAVEOUTPUT(ice1, ice2, Tsfc, Qm, Qf, Qe, Qh, Qc, chi, ...
%    balance, dt_sum, T, f_ice, f_liq, d_liq, d_lyr, d_evp, Sc, errH, errT, ...
%    iter, flag)
%
% % Calling syntax:
% % [ice1, ice2] = SAVEOUTPUT(ice1, ice2, Tsfc, Qm, Qf, Qe, Qh, Qc, ...
% %             chi, balance, dt_sum, T, f_ice, f_liq, d_liq, d_lyr, d_evp, ...
% %             Sc, errH, errT, iter, opts);
%
% % save the surface energy balance
% ice1.Tsfc(iter,1)       =  Tsfc;          % surface temp
% ice1.Qm(iter,1)         =  Qm;            % melt energy
% ice1.Qf(iter,1)         =  Qf;            % freeze deficit
% ice1.Qe(iter,1)         =  Qe;            % latent
% ice1.Qh(iter,1)         =  Qh;            % sensible
% ice1.Qc(iter,1)         =  Qc;            % conductive into surface
% ice1.chi(iter,1)        =  chi;           % skin parameter
% ice1.balance(iter,1)    =  balance;       % SEB residual
% ice1.dt(iter,1)         =  dt_sum;        % check
% %ice1.zD(iter,1)        =  zD;
%
% % Save the ice column data
% ice2.Tice(:,iter)       =  T;             % ice temperature
% ice2.f_ice(:,iter)      =  f_ice;         % fraction ice
% ice2.f_liq(:,iter)      =  f_liq;         % fraction liq
% ice2.df_liq(:,iter)     =  d_liq;
% ice2.df_lyr(:,iter)     =  d_lyr;
% ice2.df_evp(:,iter)     =  d_evp;
% ice2.Sc(:,iter)         =  Sc;            % source term
% ice2.errH(:,iter)       =  errH;          % enthalpy error
% ice2.errT(:,iter)       =  errT;          % temperature error
%
% % save the diagnostic data
% % diags.Tflag(iter,1)     =  Tflag;
% % diags.LCflag(iter,1)    =  LCflag(1);
%
%
% % Sector runs:
%
% % save the surface energy balance
% ice1.Tsfc(iter,1)    =  Tsfc;       % surface temp
%
% ice2.Tice(:,iter)    =  T;             % ice temperature
% ice2.f_ice(:,iter)   =  f_ice;         % fraction ice
% ice2.f_liq(:,iter)   =  f_liq;         % fraction liq
% ice2.df_liq(:,iter)  =  d_liq;
% ice2.df_lyr(:,iter)  =  d_lyr;
%
% % % for a stripped-down run this is all that's needed:
% % T_sfc    (iter,1)    =  Tsfc;
% % T_ice    (:,iter)    =  T;
% % frac_ice (:,iter)    =  f_ice;
% % frac_liq (:,iter)    =  f_liq;
% % df_liq   (:,iter)    =  d_liq;
