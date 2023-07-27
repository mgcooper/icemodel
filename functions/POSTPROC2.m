function [ice1,ice2] = POSTPROC2(ice1,ice2,Time,opts)

% load physical constants   
load('PHYSCONS','Tf','ro_liq','Ls','Lf');

% calculate surface and subsurface runoff
if strcmp('skinmodel', opts.simmodel)
   ice1 = SRFRUNOFF(ice1,ro_liq,Ls,Lf,opts.dt);
elseif strcmp('icemodel', opts.simmodel)
   ice1 = ICERUNOFF(ice1,ice2,opts);
end

% convert temperature to celsius
ice1.Tsfc = min(ice1.Tsfc-Tf,0);
ice2.Tice = min(ice2.Tice-Tf,0);

% convert to timetable
ice1 = struct2table(ice1);
ice1 = table2timetable(ice1,'RowTimes',Time);

% Retime to hourly
if opts.dt == 900
   
   ice1  = retime(ice1,'hourly','mean');
   feb29 = month(ice1.Time) == 2 & day(ice1.Time) == 29;
   ice1  = ice1(~feb29,:);

   % init tmp arrays to retime ice2
   tmp.Tice = nan(size(ice2.Tice,1),numel(ice1.Time));

   if strcmp('icemodel', opts.simmodel)
      tmp.df_liq = nan(size(ice2.df_liq,1),numel(ice1.Time));
      tmp.f_ice = nan(size(ice2.f_ice,1),numel(ice1.Time));
      tmp.f_liq = nan(size(ice2.f_liq,1),numel(ice1.Time));
   end

   for n = 1:numel(ice1.Time)

      % this works b/c we know it's fifteen minute data
      i1 = n*4-3;
      i2 = n*4;

      tmp.Tice(:,n) = mean(ice2.Tice(:,i1:i2),2);

      if strcmp('icemodel', opts.simmodel)
         tmp.df_liq(:,n) = sum(ice2.df_liq(:,i1:i2),2);
         tmp.f_liq(:,n) = mean(ice2.f_liq(:,i1:i2),2);
         tmp.f_ice(:,n) = mean(ice2.f_ice(:,i1:i2),2);
      end
      
   end

   ice2 = tmp;
end

% round. this should be faster than the two options below
ice1.Tsfc = round(ice1.Tsfc,5);
ice1.melt = round(ice1.melt,5);
ice1.freeze = round(ice1.freeze,5);
ice1.runoff = round(ice1.runoff,5);

ice2.Tice = round(ice2.Tice,3);

if strcmp('icemodel', opts.simmodel)
   ice2.f_ice = round(ice2.f_ice,5);
   ice2.f_liq = round(ice2.f_liq,5);
   ice2.df_liq = round(ice2.df_liq,8);
end

if isfield(ice1,'drain')
   ice1.drain = round(ice1.drain,5);
end


% Retime ice2 to hourly and round the data
% fields = fieldnames(ice2);
% oldtime = dates;
% newtime = ice1.Time;

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
