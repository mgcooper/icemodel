% function [f,ax] = plotice2(X,Y,Z,varargin)
function h = plotice2(ice2,varname,varargin)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
t1 =  datetime(year(ice2.Time(1)),7,1,0,0,0,'TimeZone','UTC');
t2 =  datetime(year(ice2.Time(1)),7,31,23,0,0,'TimeZone','UTC');

p  = MipInputParser;
p.FunctionName = 'plotice2';
p.addRequired('ice2',@(x)isstruct(x));
p.addRequired('varname',@(x)ischar(x)||isstring(x));
p.addParameter('t1',t1,@(x)isdatetime(x));
p.addParameter('t2',t2,@(x)isdatetime(x));
p.addParameter('hoursofday',1:24,@(x)isnumeric(x));
p.parseMagically('caller');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tidx     =  isbetween(ice2.Time,t1,t2);
zidx     =  1:find(roundn(ice2.Z,-1)==2,1,'last');
X        =  ice2.Time(tidx);
Y        =  ice2.Z(zidx);
V        =  ice2.(varname)(zidx,tidx);
nz       =  size(V,1);
ndays    =  size(V,2)/24;
Vrs      =  reshape(V,nz,24,ndays);
Xrs      =  reshape(X,24,ndays);
X        =  mean(Xrs(hoursofday,:),1);
V        =  squeeze(mean(Vrs(:,hoursofday,:),2));

h.figure =  figure;
h.pcolor =  pcolor(X,Y,V); shading flat;
h.ax     =  gca;
 
set(gca,'YDir','reverse','YLim',[0 2]);
ylabel('Depth [m]');
c = colorbar;
axis tight
box off

switch varname
   case 'k_eff'
      title = '$k_{eff}\;$ [W m$^{-1}$ K$^{-1}$]';
   case {'ro_snow','ro_ice'}
      title = '$\rho_{ice}\;$ [kg m$^{-3}$]';
   case 'Tice'
      title = '$T\;$ [$\circ$C]';
   case 'cp_sno'
      title = '$C_{p}\;$ [J kg$^{-1}$ K$^{-1}$]';
   case 'f_liq'
      title = '$\theta_{liq}\;$ [-]';
   case 'rel_sat'
      title = 'Relative Saturation [-]';
end

%       c.Label.String = '$k_{eff}$';
%       c.Label.Interpreter = 'latex';
   
   ct             = get(c,'title');
   ct.String      = title;
   ct.Units       = 'normalized';
   ct.Interpreter = 'latex';


% % build a grid to plot the subsurface
%    Time     =  ice2.Time;
%    ndays    =  height(Time)/24;
%    Z        =  ice2.Z;
%    V        =  ice2.(varname);
%    tidx     =  isbetween(Time,t1,t2);
%    tidx2    =  reshape(isbetween(Time,t1,t2),24,ndays);
%    zidx     =  find(roundn(Z,-1)==2,1,'last');
%    t        =  reshape(Time,24,ndays);
%    t        =  t(12,tidx2(hoursofday,:));
%    [X,Y]    =  meshgrid(t,Z(1:zidx));
%    nlayers  =  size(ice2.(var),1);


% narginchk(3,5)
% if nargin==5
%     t1 = varargin{1};
%     t2 = varargin{2};
% else
%     
% end


