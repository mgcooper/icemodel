function Tsfc = SEBSOLVE_robust(EEE,AAA,FFF,CCC,Tair,ea,wspd,emiss,SB,Tf,...
   xTsfc,scoef,fopts,liqflag,metiter)

Sfnc = @STABLEFN;
Vfnc = @VAPOR;
fSEB = @(Tsfc) EEE - emiss*SB.*Tsfc.^4 + ...
   AAA.*Sfnc(Tair,Tsfc,wspd,scoef).*(Tair-Tsfc) + ...
   FFF.*CCC.*Sfnc(Tair,Tsfc,wspd,scoef) .* ...
   (ea-Vfnc(Tsfc,Tf,liqflag));

% Note: I should add a call to SOLVE

% find the region with zero crossing
[a,b] = fzero_guess_to_bounds(fSEB,xTsfc,[xTsfc-50 xTsfc+50]);

if isnan(a)
   try
      Tsfc = fzero(fSEB,xTsfc,fopts);
   catch
      Tsfc = min(Tair,Tf);
   end
else
   try
      Tsfc = fzero_brent(fSEB, a, b, fopts.TolX);
   catch
      % fallback to fsolve
      try
         Tsfc = fzero(fSEB,[a b],fopts);
      catch 
         % fallback to secant method (or fsolve? )
         try
            Tsfc = fzero_secant(fSEB, a, b, fopts.TolX);
         catch
            Tsfc = min(Tair,Tf); % set to air temperature if all methods fail
         end
      end
   end
end

% This should only occur during spinup, and only if built-in fzero is used
if isnan(Tsfc)
   Tsfc = min(Tair,Tf);
end

%% more robust

% % find the region with zero crossing
% a = nan;
% iter = 0;
% while isnan(a) && iter<15 % stop at +/- 120K
%    iter = iter+1;
%    [a,b] = fzero_guess_to_bounds(fSEB,xTsfc,[xTsfc-(45+5*iter) xTsfc+(45+5*iter)]);
% end
%
% if isnan(a)
% %    try
% %       Tsfc = fzero(fSEB,xTsfc); % try unconstrained (not during spinup though)
% %    catch ME
%       Tsfc = Tair;
% %    end
% else
%    Tsfc = fzero_brent(fSEB,[a b],fopts.TolX);
%    % Tsfc = fzero(fSEB,[a b],fopts);
% end

% % this should only occur during spinup, and only if built-in fzero is used
% if isnan(Tsfc)
%    Tsfc = Tair;
% end
