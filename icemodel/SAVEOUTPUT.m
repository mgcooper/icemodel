function [ice1, ice2] = SAVEOUTPUT(step, ice1, ice2, vars1, vars2, data1, data2)
   %SAVEOUTPUT
   %
   %#codegen

   % Note: the order of vars1 and vars2 (set in icemodel.setopts) must match the
   % order of data1, data2 (passed in from icemodel main). It's fine if more
   % variables exist in vars1, vars2 than there are data in data1, data2, but
   % the vars1, vars2 lists must be in the same order as the data in data1,
   % data2. Thus vars1, vars2 can be set to some base set of variables in
   % setopts, and the actual variables can be adjusted in icemodel main.
   for n = 1:numel(data1)
      ice1.(vars1{n})(step, 1) = data1{n};
   end

   for n = 1:numel(data2)
      ice2.(vars2{n})(:, step) = data2{n};
   end
end

% function [ice1, ice2, diags] = SAVEOUTPUT(step, ice1, ice2, vars1, vars2, ...
%       data1, data2, diags1, diags2, diagvars1, diagvars2, diagdata1, diagdata2)
%
%    if nargin > 7
%       for n = 1:numel(diagdata1)
%          diags1.(diagvars1{n})(step, 1) = diagdata1{n};
%       end
%       for n = 1:numel(diagdata2)
%          diags2.(diagvars2{n})(:, step) = diagdata2{n};
%       end
%    end
% end

% % or:
% function [ice1, ice2] = SAVEOUTPUT(step, ice1, ice2, vars1, vars2, varargin)
%
% data1 = varargin(1:numel(vars1));
% data2 = varargin(numel(vars1)+1:numel(varargin));
%
% for n = 1:numel(vars1)
%    ice1.(vars1{n})(step, 1) = data1{n};
% end
% for n = 1:numel(vars2)
%    ice2.(vars2{n})(:, step) = data2{n};
% end
%
%
% function [ice1, ice2] = SAVEOUTPUT(ice1, ice2, Tsfc, Qm, Qf, Qe, Qh, Qc, chi, ...
%    balance, dt_sum, T, f_ice, f_liq, d_liq, d_lyr, d_evp, Sc, errH, errT, ...
%    step, flag)
%
% % Calling syntax:
% % [ice1, ice2] = SAVEOUTPUT(ice1, ice2, Tsfc, Qm, Qf, Qe, Qh, Qc, ...
% %             chi, balance, dt_sum, T, f_ice, f_liq, d_liq, d_lyr, d_evp, ...
% %             Sc, errH, errT, step, opts);
%
% % save the surface energy balance
% ice1.Tsfc(step,1)       =  Tsfc;          % surface temp
% ice1.Qm(step,1)         =  Qm;            % melt energy
% ice1.Qf(step,1)         =  Qf;            % freeze deficit
% ice1.Qe(step,1)         =  Qe;            % latent
% ice1.Qh(step,1)         =  Qh;            % sensible
% ice1.Qc(step,1)         =  Qc;            % conductive into surface
% ice1.chi(step,1)        =  chi;           % skin parameter
% ice1.balance(step,1)    =  balance;       % SEB residual
% ice1.dt(step,1)         =  dt_sum;        % check
% %ice1.zD(step,1)        =  zD;
%
% % Save the ice column data
% ice2.Tice(:,step)       =  T;             % ice temperature
% ice2.f_ice(:,step)      =  f_ice;         % fraction ice
% ice2.f_liq(:,step)      =  f_liq;         % fraction liq
% ice2.df_liq(:,step)     =  d_liq;
% ice2.df_lyr(:,step)     =  d_lyr;
% ice2.df_evp(:,step)     =  d_evp;
% ice2.Sc(:,step)         =  Sc;            % source term
% ice2.errH(:,step)       =  errH;          % enthalpy error
% ice2.errT(:,step)       =  errT;          % temperature error
%
% % save the diagnostic data
% % diags.Tflag(step,1)     =  Tflag;
% % diags.LCflag(step,1)    =  LCflag(1);
%
%
% % Sector runs:
%
% % save the surface energy balance
% ice1.Tsfc(step,1)    =  Tsfc;       % surface temp
%
% ice2.Tice(:,step)    =  T;             % ice temperature
% ice2.f_ice(:,step)   =  f_ice;         % fraction ice
% ice2.f_liq(:,step)   =  f_liq;         % fraction liq
% ice2.df_liq(:,step)  =  d_liq;
% ice2.df_lyr(:,step)  =  d_lyr;
%
% % % for a stripped-down run this is all that's needed:
% % T_sfc    (step,1)    =  Tsfc;
% % T_ice    (:,step)    =  T;
% % frac_ice (:,step)    =  f_ice;
% % frac_liq (:,step)    =  f_liq;
% % df_liq   (:,step)    =  d_liq;
