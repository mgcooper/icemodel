%--------------------------------------------------------------------------
%   Calculate solar heating within the ice matrix
%--------------------------------------------------------------------------

function [Sc,Sp] = SOLARHEAT(Qsi,total_solar,dQp,dz,chi)
%--------------------------------------------------------------------------

% Get the solar radiation penetrating the ice surface.
    Qsip        =   SOLARPEN(Qsi,chi);

% Scale the solar radiation absorbed in each c.v. (dQp), which is given in
%   terms of total_solar, to the observed solar radiation 
    dQp         =   (-Qsip/total_solar).*dQp;

% Compute the source term in each c.v. at the current timestep (dQ/dz)  
    Sc          =   dQp./dz;              % [W/m3]
    Sp          =   0.0.*Sc;

%     figure; loglog(dQp,0.01:0.01:12); %semilogx(dQp,0.01:0.01:12); 
%     set(gca,'YDir','reverse','YLim',[0 4])
%     xlabel('Solar radiation'); ylabel('depth (m)');
    
    