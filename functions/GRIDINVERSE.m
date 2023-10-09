%--------------------------------------------------------------------------
%   Transform absorbed solar from the spectral grid to the thermal grid
%--------------------------------------------------------------------------

function dQp   =  GRIDINVERSE(xynet,dz_spect,dz_therm,JJ_spect,JJ_therm)

% extrapolate xynet to the bottom of the thermal grid
   JJnew             =  JJ_therm*dz_therm/dz_spect;
   JJext             =  round(JJnew-JJ_spect, 0);
   dxy               =  xynet(JJ_spect)-xynet(JJ_spect-1);
   xyextrap          =  zeros(JJext,1);
   xyextrap(1)       =  xynet(JJ_spect) + dxy;
   xyextrap(2:JJext) =  xyextrap(1:JJext-1) + dxy;
   xynew             =  [xynet;xyextrap];
% aggregate it to the scale of the thermal grid
   dq                =  xynew(1:JJnew-1)-xynew(2:JJnew);
   dq(JJnew)         =  dq(JJnew-1);
% get the number of grid cells in the first 3 m segment on each grid
   rz                =  JJnew/JJ_therm;
% reshape the radiation into equal chunks along the thermal grid
   dq                =  reshape(dq,rz,JJ_therm);    
% sum up those equal chunks to get the absorbed radiation in each thermal c.v.
   dQp               =  transpose(sum(dq,1)); % == sumdQdz

%--------------------------------------------------------------------------
%   A note on what's happening here. The two-stream solution gives the
%   total upflux (X) and downflux (Y) from the respective boundary to each
%   layer interface. xynet is energy conservation applied to each c.v., in
%   the following manner: 
%
%   Y↓(1)  X↑(1)    Xnet = X↑(2)-X↑(1), Ynet = Y↓(1)-Y↓(2)
%   ____________
%   ____________    XYnet = Xnet - Ynet = (X↑(2)-X↑(1))-(Y↓(1)-Y↓(2))
%   Y↓(2)  X↑(2)  

%   So XYnet is the net flux and dq = d/dz(XYnet) is the absorbed flux. In
%   this function I first extrapolate xynet onto the thermal grid in case
%   they don't have the same extent, and then I take dq = d/dz(XYnet) to
%   get the absorbed flux. Note the way this is implemented in the code is
%   somewhat different than the diagram above, since xynet is computed in
%   SOLVETWOSTREAM2 by converting the interface fluxes to layer averaged
%   values (note the /2.0 averaging) and then the differencing is applied
%   in here to get the interface absorption, Qz. Also note that dq is scaled by
%   total_solar, and should be thought of as a percent absorption in each
%   layer. In the thermal model, SOLAR_HEAT converts dq to the actual
%   amount of absorbed solar radiation in each layer using the conversion
%   formula dQ = (-Qsip/total_solar).*dQ
%--------------------------------------------------------------------------
      
