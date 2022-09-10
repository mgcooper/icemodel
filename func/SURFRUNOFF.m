%--------------------------------------------------------------------------
%   Post process the surface runoff to account for refreezing
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function runoff = SURFRUNOFF(ice1,opts)
%--------------------------------------------------------------------------

% this says if any melt occured during the prior tlag hours then subtract
% the freeze from the melt + cond to get runoff, otherwise let melt
% note: need to add rain

tlag                =   opts.tlagsurf;

melt                =   ice1.surf_melt;
freeze              =   ice1.surf_freeze;
cond                =   ice1.surf_cond;

% to get this in the model, replace meltdiffs etc. with imelt, etc.
meltdiffs           =   [0.0 ; diff(melt)];
freezediffs         =   [0.0 ; diff(freeze)];
conddiffs           =   [0.0 ; diff(cond)];

if tlag == 0
    runoff(1,1)         =   0.0;
else
    runoff(1:tlag,1)    =   0.0;
end


for n = 1+tlag:length(melt)
    meltsumlag      =   sum(meltdiffs(n-tlag:n) + conddiffs(n-tlag:n));
    potrunoff       =   conddiffs(n) + meltdiffs(n);
    potfreeze       =   min(freezediffs(n),meltsumlag);
    potfreeze       =   max(potfreeze,0.0);
    if opts.skinfreeze == true
        if meltsumlag > 0.0
            netrunoff   =   potrunoff-potfreeze;
            runoff(n,1) =   max(runoff(n-1)+netrunoff,0.0);
        else
            runoff(n,1) =   runoff(n-1) + potrunoff;
        end 
    else
        runoff(n,1) =   runoff(n-1) + potrunoff;
    end
end

