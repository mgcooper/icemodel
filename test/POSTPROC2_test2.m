function [Tsfc,runoff,melt,freeze,                                      ...
         T_ice,f_ice,f_liq,df_liq] =   POSTPROC2_test2(T_sfc,T_ice,     ...
                                       frac_ice,frac_liq,df_liq,met,opts)

% v2 keeps everything as arrays
dz       = opts.dz_thermal;
T_sfc    = min(T_sfc-273.15,0);
T_ice    = min(T_ice-273.15,0);

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

melt        = tocolumn(cumsum(melt));     % cumulative melt
freeze      = tocolumn(cumsum(freeze));   % cumulative freeze

% Retime met to hourly

% Retime ice1 to hourly
dates          = met.Time; % dates.TimeZone = 'UTC';
% ice1           = struct2table(ice1);
% ice1           = table2timetable(ice1,'RowTimes',dates);
% ice1           = retime(ice1,'hourly','mean');

% Retime ice2 to hourly and round the data
oldtime     = dates;
newtime     = dates(1):hours(1):dates(end)-minutes(45);

% round and retime ice2 vars
T_ice    = round(transpose(interp1(oldtime,transpose(T_ice),newtime)),3);
f_ice    = round(transpose(interp1(oldtime,transpose(frac_ice),newtime)),5);
f_liq    = round(transpose(interp1(oldtime,transpose(frac_liq),newtime)),5);
df_liq   = round(transpose(interp1(oldtime,transpose(df_liq),newtime)),8);

% round and retime ice1 vars
Tsfc        = round(interp1(oldtime,T_sfc,newtime),5);
melt        = round(interp1(oldtime,melt,newtime),5);
freeze      = round(interp1(oldtime,freeze,newtime),5);
runoff      = round(interp1(oldtime,runoff,newtime),5);
