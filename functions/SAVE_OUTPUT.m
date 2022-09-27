% this isn't implemented, but might speed-up the code

function [  enbal,                                              ...
            ice1,                                               ...
            ice2    ]   =   SAVE_OUTPUT(n,iyear,imonth,iday,    ...
                            ihour,Tair,Tsfc,rh,windspd,Qsi,     ...
                            Qli,Qle,Qh,Qe,Qc,Qm,balance,        ...
                            albedo,qsfactor,Tflag,depth_melt,   ...
                            depth_freeze,depth_drain,surf_melt, ...
                            ice_surf,gamma,Cp_snow,ro_snow,     ...
                            ro_snow_dry,h_melt,h_freeze,h_drain,...
                            h_bal,Tf,T_old,enbal,ice1,ice2)
                        
      % save the desired quantities
        Tair            =   Tair-Tf;
        Tsfc            =   Tsfc-Tf;
        Tice            =   T_old-Tf;
        Qsip            =   (1-qsfactor)*(1-albedo)*Qsi;
        Qsrf            =   (qsfactor)*(1-albedo)*Qsi;
        ro_snow1        =   ro_snow_dry(1);
        Tice1           =   T_old(1)-Tf;

      % save surface energy balance
        enbal_n         =   table(n,iyear,imonth,iday,ihour,            ...
                            Tair,Tsfc,rh,windspd,Qsi,Qli,Qle,           ...
                            Qh,Qe,Qc,Qm,Qsip,Qsrf,balance,              ...
                            albedo,Tflag);
        
        ice1_n          =   table(n,iyear,imonth,iday,ihour,            ...
                            Tair,Tsfc,Tice1,depth_melt,                 ...
                            depth_freeze,depth_drain,surf_melt,         ...
                            ice_surf,ro_snow1);
                        
      % concatenate today's data to the tables
        enbal           =   [enbal;enbal_n];
        ice1            =   [ice1;ice1_n];

        % use a structure for the vertical column data
        ice2.Tice(:,n)      =   Tice;
        ice2.gamma(:,n)     =   gamma;
        ice2.Cp_snow(:,n)   =   Cp_snow;
        ice2.ro_snow(:,n)   =   ro_snow;
        ice2.ro_snow_d(:,n) =   ro_snow_dry;
        ice2.h_melt(:,n)    =   h_melt;
        ice2.h_freeze(:,n)  =   h_freeze;
        ice2.h_drain(:,n)   =   h_drain;
        ice2.h_bal(:,n)     =   h_bal;

end
    

