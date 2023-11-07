function [f_ice, f_liq, d_drn, xsubl] = ICESUBL(f_ice, f_liq, d_drn, ro_ice, ...
      ro_liq, Ls, Lv, f_ice_min, f_liq_min, dz_therm, dt_new, Qe, liqflag)
   %SUBL Compute ice sublimation from top layer control volume
   %
   % f_ice = fraction of ice by volume in each control volume
   % f_liq = fraction of liquid water by volume in each control volume
   % d_drn = change in f_liq due to water drainage from control volume
   % xsubl = "extra" sublimation i.e., error due to sublimation exceeding f_ice
   %
   %#codegen

   % not sure about liqflag and f_liq_min because liqflag means f_liq>0.02, so
   % once f_liq gets below that, no evaporation will occur, also, condensation
   % should be able to occur when liqflag is false

   % potential evaporation
   pevap = Qe / (Lv * ro_liq) * dt_new / dz_therm;

   % top-layer liquid water and ice fraction
   f_liq_top = f_liq(1);
   f_ice_top = f_ice(1);

   % Initialize variables
   xevap = 0;
   xsubl = 0;

   % note: if liqflag is true, f_liq must be > f_liq_min, but if liqflag is
   % false and pevap > 0, it is possible that f_liq is < f_liq_min on entry to
   % this function. CORRECTION - liqflag can be true and f_liq can be <
   % f_liq_min on entry to this function b/c f_liq_min is based on the start of
   % the timestep f_wat, and during ICEENBAL, f_liq(1) can drop below f_liq_min
   % as long as these evaluate false:
   % f_liq_top < 0.95*f_liq_min % || f_liq_top > 1.05*f_liq_max

   if liqflag == true || pevap > 0

      if pevap < 0 % evaporation

         if f_liq_top < f_liq_min % f_liq_top - f_liq_min < 0
            % only residual water exists, send pevap to SUBL
            xevap = pevap;

            % add metiter input to use this but it never triggered in tests
            %fprintf('metiter = %d, f_liq(1) < f_liq_min\n',metiter)

         elseif abs(pevap) <= (f_liq_top - f_liq_min)
            % some water evaporates
            % d_evap = pevap;
            f_liq(1) = f_liq(1) - abs(pevap); % aevap = pevap, xevap = 0

         else
            % all available water evaporates
            aevap = -(f_liq_top - f_liq_min);
            f_liq(1) = f_liq_min; % keep residual water
            xevap = pevap - aevap; % extra evap becomes potential subl
         end

         % send demand not satisfied by evaporation to SUBL
         [f_ice, xsubl] = SUBL(xevap, f_ice, f_ice_min, ro_ice, ro_liq, Lv, Ls);

      elseif pevap > 0 % condensation

         if pevap <= (1 - f_liq_top - f_ice_top)
            % all condensation stored
            f_liq(1) = f_liq(1) + pevap;  % aevap = pevap, xevap = 0

         elseif pevap > (1 - f_liq_top - f_ice_top)
            % some condensation stored, some converts to runoff
            aevap = 1 - f_liq_top - f_ice_top;
            f_liq(1) = 1 - f_ice_top;
            d_drn(1) = d_drn(1) + (pevap - aevap); % xevap = pevap - aevap
         end
      end

   else % sublimation
      [f_ice, xsubl] = SUBL(pevap, f_ice, f_ice_min, ro_ice, ro_liq, Lv, Ls);
   end
end

%% Sub functions
function [f_ice, xsubl] = SUBL(pevap, f_ice, f_ice_min, ro_ice, ro_liq, Lv, Ls)

   xsubl = 0;
   if pevap >= 0
      return
   end

   % convert evap to ice frac-equivalent
   psubl = pevap * (ro_liq * Lv) / (ro_ice * Ls);
   f_ice_top = f_ice(1);

   % Note: since layer combination is based on f_ice, I don't think we need to
   % check against f_min, we only need the < 0 check. If f_ice(1)+psubl is < 0,
   % then we have an error, otherwise if f_ice(1)+psubl < f_min, the layers will
   % be combined.

   if f_ice_top < f_ice_min

      if abs(psubl) < f_ice_top
         % some ice sublimates, and the top layer will be combined in ICEMF
         f_ice(1) = f_ice(1) - abs(psubl);

      elseif abs(psubl) >= f_ice_top
         % all ice sublimates (asubl = f_ice(1) )
         f_ice(1) = 0.0;

         % "extra" sublimation that cannot be satisfied
         xsubl = psubl + f_ice_top;
      end

   elseif abs(psubl) < (f_ice_top - f_ice_min)
      % some ice sublimates (asubl = psubl)
      f_ice(1) = f_ice(1) - abs(psubl);

   elseif abs(psubl) >= (f_ice_top - f_ice_min)
      % all ice sublimates (asubl = -(fi - f_ice_min) )
      f_ice(1) = f_ice_min;         % keep minimum ice thickness

      % "extra" sublimation that cannot be satisfied
      xsubl = psubl + f_ice_top - f_ice_min;
   end
end
% % might add this:
% subl = @(e) e .* roLv./roLs ;

% but in icemodel, roL is set to roLv or roLs depending on liqflag, so Qe is
% already computed wrt to them. That way evap/subl are computed using the same
% formula: e = Qe/roL*dt/dz. The conversion is only needed when there is "extra
% evap", i.e., if evap removes all f_liq above f_min, and the remainder of Qe is
% allocated to subl, then we use the conversion: s = xe*roLv/roLs;

%%

% if liqflag == true || pevap > 0
%
%    if pevap < 0
%
%       if f_liq_top < f_liq_min

% % Notes here are for the case above where liqflag is true, but once we get to
% ICEMF, f_liq(1) is < f_liq_min due to refreezing within the timestep, below is
% an option to handle this that respects liqflag true, but instead I think it's
% better to instead do sublimation, and convert the liqflag true latent heat
% from vaporization to sublimation

% fprintf('metiter = %d, f_liq(1) < f_liq_min\n',metiter)
% % There are two options to deal with this edge case.
% % 1) The routine below which should work in all cases. If "some water
% % evaporates", then f_liq will drop further below f_liq_min, and on
% % return to ICEMF, T(1) will be updated to be consistent with the new
% % f_liq(1), which should prevent triggering an error in MZTRANSFORM on
% % the next timestep ... but pay attention to the possibility that if
% % f_liq(1) goes below f_liq_min, depending on the T(1) update, it might
% % end up triggering in an infinite recursion in ICEENBAL b/c restarting
% % the timestep won't work no matter how low we go on dt_new, but the
% % more I think about it, it should be fine, so pretend it works, then
% % the important thing is that the call to SUBL(xevap, ...) must be
% % moved out of the else block, and into a new "if xevap>0" block at the
% % end (or just call SUBL with xevap = 0 and it will return).
% %
% % 2) Call SUBL instead
%
% if abs(pevap) < f_liq_top
%    % some water evaporates
%    f_liq(1) = f_liq(1) - abs(pevap);
%
% elseif abs(pevap) >= f_liq_top
%    % all water evaporates (aevap = f_liq(1) )
%    f_liq(1) = 0.0;
%
%    % "extra" evaporation that cannot be satisfied
%    xevap = pevap + f_liq_top;
% end
%
% % for debugging
% % 273.150-sqrt(((0.0007895164+f_ice_top.*0.9170)./0.0007895164-1))./100
