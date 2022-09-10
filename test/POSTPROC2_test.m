function [ice1,ice2] = POSTPROC2_test(T_sfc,T_ice,frac_ice,frac_liq,df_liq,met,opts)

% this one keeps ice1 as a struct. v2 keeps everything as arrays

dz             = opts.dz_thermal;
ice1.Tsfc      = min(T_sfc-273.15,0);
ice2.df_liq    = df_liq;
ice2.f_ice     = frac_ice;
ice2.f_liq     = frac_liq;
ice2.Tice      = min(T_ice-273.15,0);

% partition runoff into melt/freeze
melt     = zeros(size(T_sfc));
freeze   = zeros(size(T_sfc));

for n = 1:numel(T_sfc)
    imlt        = df_liq(:,n)>0;
    ifrz        = df_liq(:,n)<0;
    melt(n)     = sum(dz(1).*df_liq(imlt,n));
    freeze(n)   = sum(-dz(1).*df_liq(ifrz,n));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute runoff using a time lag

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

ice1.runoff       = runoff;                     % cumulative runoff
ice1.melt         = tocolumn(cumsum(melt));     % cumulative melt
ice1.freeze       = tocolumn(cumsum(freeze));   % cumulative freeze

% Retime ice1 to hourly
dates          = met.Time;

% Retime ice2 to hourly and round the data
fields1        = fieldnames(ice1);
fields2        = fieldnames(ice2);
oldtime        = dates;
newtime        = dates(1):hours(1):dates(end)-minutes(45);

for mm = 1:numel(fields1)

   thisfield         = fields1{mm};
   newdat            = interp1(oldtime,ice1.(thisfield),newtime);
   ice1.(thisfield)  = transpose(newdat);
   ice1.(thisfield)  = round(ice1.(thisfield),5);

end

for mm = 1:numel(fields2)

   thisfield         = fields2{mm};
   newdat            = interp1(oldtime,transpose(ice2.(thisfield)),newtime);
   ice2.(thisfield)  = transpose(newdat);

   if ismember(thisfield,{'f_ice','f_liq'})
      ice2.(thisfield) = round(ice2.(thisfield),5);
   elseif ismember(thisfield,{'Tice'})
      ice2.(thisfield) = round(ice2.(thisfield),3);
   elseif ismember(thisfield,{'df_liq'})
      ice2.(thisfield) = round(ice2.(thisfield),8);
   end
end