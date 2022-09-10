function [ice1,ice2] = POSTPROC2(T_sfc,T_ice,frac_ice,frac_liq,df_liq,Time,opts)
   
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

% compute runoff using a time lag
runoff         = COLUMNRUNOFF(melt,freeze,opts);

ice1.runoff    = runoff;                     % cumulative runoff
ice1.melt      = tocolumn(cumsum(melt));     % cumulative melt
ice1.freeze    = tocolumn(cumsum(freeze));   % cumulative freeze

% Retime ice1 to hourly and round
ice1           = struct2table(ice1);
ice1           = table2timetable(ice1,'RowTimes',Time);
ice1           = retime(ice1,'hourly','mean');
ice1           = rmleapinds(ice1);

% init temp arrays
tmp.df_liq     = nan(size(ice2.df_liq,1),numel(ice1.Time));
tmp.f_ice      = nan(size(ice2.f_ice,1),numel(ice1.Time));
tmp.f_liq      = nan(size(ice2.f_liq,1),numel(ice1.Time));
tmp.Tice       = nan(size(ice2.Tice,1),numel(ice1.Time));

for n = 1:numel(ice1.Time)

   % this works b/c we know it's fifteen minute data
   i1    = n*4-3;
   i2    = n*4;
   
   tmp.df_liq(:,n)   =  sum(ice2.df_liq(:,i1:i2),2);
   tmp.f_liq(:,n)    =  mean(ice2.f_liq(:,i1:i2),2);
   tmp.f_ice(:,n)    =  mean(ice2.f_ice(:,i1:i2),2);
   tmp.Tice(:,n)     =  mean(ice2.Tice(:,i1:i2),2);
end

ice2 = tmp;

% this should be faster than the two options below
ice1.Tsfc      = round(ice1.Tsfc,5);
ice1.runoff    = round(ice1.runoff,5);
ice1.melt      = round(ice1.melt,5);
ice1.freeze    = round(ice1.freeze,5);
ice2.df_liq    = round(ice2.df_liq,8);
ice2.f_ice     = round(ice2.f_ice,5);
ice2.f_liq     = round(ice2.f_liq,5);
ice2.Tice      = round(ice2.Tice,3);


% Retime ice2 to hourly and round the data
% fields         = fieldnames(ice2);
% oldtime        = dates;
% newtime        = ice1.Time;

% for mm = 1:numel(fields)
% 
%    thisfield         = fields{mm};
%    newdat            = interp1(oldtime,transpose(ice2.(thisfield)),newtime);
%    ice2.(thisfield)  = transpose(newdat);
% 
%    if ismember(thisfield,{'f_ice','f_liq'})
%       ice2.(thisfield) = round(ice2.(thisfield),5);
%    elseif ismember(thisfield,{'Tice'})
%       ice2.(thisfield) = round(ice2.(thisfield),3);
%    elseif ismember(thisfield,{'df_liq'})
%       ice2.(thisfield) = round(ice2.(thisfield),8);
%    end
% end

% repeat rounding for ice1. in this case, we round all fields to 5
% digits, so i use 'fields' in the logical checks below, but if i had
% vars i wanted rounded to different precision, would need to replace 

% fields = ice1.Properties.VariableNames;
% for mm = 1:numel(fields)
%    thisfield = fields{mm};      
%    if ismember(thisfield,fields)
%       ice1.(thisfield) = round(ice1.(thisfield),5);
%    end
% end

% ice2.Z         = (dz/2:dz:Z-dz/2)';
