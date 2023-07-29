function [  f_liq_j,                                                    ...
            f_ice_j,                                                    ...
            T_j,                                                        ...
            Sc_j,                                                       ...
            Sp_j ]   =  COMBINEHEAT(f_liq,f_ice,T,Sc,Sp,Tf,TL,          ...
                        cp_ice,Lf,ro_ice,ro_liq,fcp,j1,j2,fopts)
%COMBINEHEAT combine layers, conserving mass, enthalpy, and absorbed radiation

% j1 is the thin layer that disappears. After the combination is finished,
% j1 is always removed and the combined values get put into j2

% pay attention to g_watC vs g_wat1+g_wat2, they are related as:
% g_watC = (g_wat1 + g_wat2)/2, so long as dz is constant.

% combined solar heat
Sc_j = Sc(j1) + Sc(j2);
Sp_j = Sp(j1) + Sp(j2);

% combined mass
g_liq12 = ro_liq*(f_liq(j1)+f_liq(j2));
g_wat1 = f_liq(j1)*ro_liq + f_ice(j1)*ro_ice;      % [kg m-3]
g_wat2 = f_liq(j2)*ro_liq + f_ice(j2)*ro_ice;
f_wat1 = f_liq(j1) + f_ice(j1)*ro_ice/ro_liq;      % frac_wat_old
f_wat2 = f_liq(j2) + f_ice(j2)*ro_ice/ro_liq;      % frac_wat_old
f_watC = (f_wat1 + f_wat2)/2;

% depression temp's
Td1 = Tf-T(j1);
Td2 = Tf-T(j2);

% combined temperature of dry and wet snow
if T(j1) == T(j2)
   T_j = T(j1);                                 % case 1: equal
   TdC = Td1;
elseif T(j1) < TL && T(j2) < TL                    % case 2: !equal, dry
   T_j = (T(j1)*g_wat1 + T(j2)*g_wat2)/(g_wat1+g_wat2);
   TdC = Tf-T_j;
else                                               % case 3: !equal, wet
   A = fcp^2*cp_ice*(g_wat1+g_wat2);
   B = fcp^2*(Lf*g_liq12-cp_ice*(g_wat1*Td1+g_wat2*Td2));
   C = A/fcp^2;
   D = B/fcp^2 - Lf*(g_wat1+g_wat2);
   %TdC = roots([A B C D]);       % should be second root
   fTdC = @(TdC) A*TdC^3 + B*TdC^2 + C*TdC + D;

   % % to make A = 1:
   % A = 1;
   % B = (Lf*g_liq12-cp_ice*(g_wat1*Td1+g_wat2*Td2))/(cp_ice*(g_wat1+g_wat2));
   % C = 1/fcp^2;
   % D = (Lf*(g_liq12-(g_wat1+g_wat2))-cp_ice*(g_wat1*Td1+g_wat2*Td2))/...
   %          (fcp^2*cp_ice*(g_wat1+g_wat2));

   % try an endpoint constrained range
   try
      % [a,b] = fzero_guess_to_bounds(fTdC,min(Td1,Td2),max(Td1,Td2));
      % TdC = fzero_brent(fTdC,a,b,fopts.TolX);
      TdC = fzero(fTdC,[min(Td1,Td2) max(Td1,Td2)],fopts);

   catch e %#ok<*NASGU>

      if strcmp(e.identifier,'MATLAB:fzero:ValuesAtEndPtsSameSign') || ...
            strcmp(e.identifier,'MATLAB:fzero:Arg2NotFinite')
         m = 'Layer combination failed using endpoints, using midpoint instead';
         e = addCause(e, MException('icemodel:COMBINEHEAT:rootfinding', m));
      end

      % try a midpoint
      try
         TdC = fzero(fTdC,(Td1+Td2)/2,fopts);
      catch e
         if strcmp(e.identifier,'MATLAB:fzero:ValuesAtEndPtsSameSign') || ...
               strcmp(e.identifier,'MATLAB:fzero:Arg2NotFinite')
            m = 'Layer combination failed using midpoint';
            e = addCause(e, MException('icemodel:COMBINEHEAT:rootfinding', m));
         end

         % use the average temperature
         try
            T_j = (T(j1)*g_wat1 + T(j2)*g_wat2)/(g_wat1+g_wat2);
            TdC = Tf-T_j;
         catch
            % If it fails here, throw the error. If this ever happens, an
            % alternative is to invoke the cv-balance subroutine to expel
            % liquid water and combine the dry ice masses, but in testing 
            % this has never happened.
            rethrow(e); 
         end

      end
   end
   T_j = Tf-TdC;
end

% If case 3 produces a fully-melted node with T>Tf, this will set T=Tf
% and ICEMF will remove the node
if TdC < 0
   T_j = (T(j1)*g_wat1 + T(j2)*g_wat2)/(g_wat1+g_wat2);
   TdC = Tf-T_j;
end

% Invert T_j to get f_liq
f_liq_j = f_watC./(1.0+(fcp.*TdC).^2.0);      % eq 67, Jordan
f_ice_j = (f_watC-f_liq_j)*ro_liq/ro_ice;
