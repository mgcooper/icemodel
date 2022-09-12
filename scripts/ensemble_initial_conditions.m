clean

savedata    =   true;
vers        =  'v10';
run         =  'run1';
pathdata    =  ['/Users/coop558/mydata/icemodel/output/' vers '/' run '/'];
pathsave    =  setpath(['GREENLAND/icemodel/input/init/' vers '/' run '/']);
patheval    =  '/Users/coop558/mydata/icemodel/eval/v11/run2/';


load([pathdata 'icemodel_ak4_2015_mar_swap_NONE_albedo.mat'],'ice2');
f_ice1   = ice2.f_ice(:,1);
f_liq1   = ice2.f_liq(:,1);    
T1       = ice2.Tice(:,1) + 273.15;
cp_sno1  = ice2.cp_sno(:,1);
k_eff1   = ice2.k_eff(:,1);
ro_sno1  = ice2.ro_sno(:,1);

load([pathdata 'icemodel_behar_2015_mar_swap_NONE_albedo.mat'],'ice2');
f_ice2   = ice2.f_ice(:,1);
f_liq2   = ice2.f_liq(:,1);    
T2       = ice2.Tice(:,1) + 273.15;
cp_sno2  = ice2.cp_sno(:,1);
k_eff2   = ice2.k_eff(:,1);
ro_sno2  = ice2.ro_sno(:,1);

% average them
init.f_ice    = (f_ice1+f_ice2)./2;
init.f_liq    = (f_liq1+f_liq2)./2;
init.T        = (T1+T2)/2;
init.cp_sno   = (cp_sno1+cp_sno2)./2;
init.k_eff    = (k_eff1+k_eff2)./2;
init.ro_sno   = (ro_sno1+ro_sno2)./2;

if savedata == true
   if ~exist(pathsave,'dir'); mkdir(pathsave); end
   save([pathsave 'init'],'init');
end

