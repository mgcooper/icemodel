%--------------------------------------------------------------------------
%   Post process the surface runoff to account for refreezing
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function runoff = COLUMNRUNOFF(melt,freeze,opts)
%--------------------------------------------------------------------------

tlag        =  opts.tlagcolumn;
runoff      =  zeros(size(melt));

for n = 1+tlag:length(melt)
    meltsumlag      =   sum(melt(n-tlag:n));
    potrunoff       =   melt(n);
    potfreeze       =   min(freeze(n),meltsumlag);
    potfreeze       =   max(potfreeze,0.0);
    if meltsumlag > 0.0
        netrunoff   =   potrunoff-potfreeze;
        runoff(n,1) =   max(runoff(n-1)+netrunoff,0.0);
    else
        runoff(n,1) =   runoff(n-1) + potrunoff;
    end 
end

% This function subtracts freeze from melt to get runoff.
% meltsumlag is the sum of melt production during the prior tlag hours up
% to and including the present hour, so tlag = 1 means you sum up 2 hours
% worth of melt, tlag = 0 means you just get the current hour of melt.
% potential freeze (potfreeze) cannot exceed this amount. In this way you
% can make an assumption about the residence time of liquid water in your
% system, and prevent freezing of water outside of that timescale. Note
% that the potential freeze and potential runoff are for the current
% timestep, so once there is no freeze energy, no water can be frozen (i.e.
% freeze energy from the middle of the night cannot be 'carried forward' or
% anything like that, but melt production from prior timesteps can be
% carried forward and subjected to the freeze energy at later timesteps,
% just like it occurs in the real world). 
% note: need to add rain