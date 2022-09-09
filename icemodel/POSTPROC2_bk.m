function [ice1,ice2] = POSTPROC2(T_sfc,T_ice,frac_ice,frac_liq,df_liq,met,opts)
   
dz             = opts.dz_thermal;
ice1.Tsfc      = min(T_sfc-273.15,0);
ice2.df_liq    = df_liq;
ice2.f_ice     = frac_ice;
ice2.f_liq     = frac_liq;
ice2.Tice      = min(T_ice-273.15,0);
      
      
runoff   = tocolumn(cumsum( sum( dz.*df_liq ) ));
   
% partition runoff into melt/freeze
melt     = zeros(size(runoff));
freeze   = zeros(size(runoff));
h_melt   = zeros(size(df_liq));
h_freeze = zeros(size(df_liq));

for n = 1:numel(runoff)
    imlt        = df_liq(:,n)>0;
    ifrz        = df_liq(:,n)<0;
    melt(n)     = sum(dz(1).*df_liq(imlt,n));
    freeze(n)   = sum(-dz(1).*df_liq(ifrz,n));
    
    % this mimics the h_freeze and h_melt that used to be in ice2
    h_melt(imlt,n)     = dz.*df_liq(imlt,n);
    h_freeze(ifrz,n)   = -dz.*df_liq(ifrz,n);
end

ice1.runoff          = runoff;                     % cumulative runoff
ice1.depth_melt      = tocolumn(cumsum(melt));     % cumulative melt
ice1.depth_freeze    = tocolumn(cumsum(freeze));   % cumulative freeze
ice1.surf_runoff     = ice1.depth_melt;
ice1.column_runoff   = COLUMNRUNOFF(h_melt,h_freeze,opts);

% Retime ice1 to hourly
dates          = met.Time; dates.TimeZone = 'UTC';
ice1           = struct2table(ice1);
ice1           = table2timetable(ice1,'RowTimes',dates);
ice1           = retime(ice1,'hourly','mean');

% Retime ice2 to hourly and round the data
fields         = fieldnames(ice2);
oldtime        = dates;
newtime        = ice1.Time;

for mm = 1:numel(fields)

   thisfield         = fields{mm};
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

% repeat rounding for ice1. in this case, we round all fields to 5
% digits, so i use 'fields' in the logical checks below, but if i had
% vars i wanted rounded to different precision, would need to replace 
fields = ice1.Properties.VariableNames;
for mm = 1:numel(fields)
   thisfield = fields{mm};      
   if ismember(thisfield,fields)
      ice1.(thisfield) = round(ice1.(thisfield),5);
   end
end

% ice2.Z         = (dz/2:dz:Z-dz/2)';
