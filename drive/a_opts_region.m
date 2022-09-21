function opts = a_opts_region(opts,sitename,meltmodel,startyear,endyear)
   
   opts.sitename        =  sitename;
   
   % model options
   opts.annual_loops    =  1;
   opts.calendar_type   =  'noleap';
   opts.dt              =  900.0;  % timestep                          [s]
   opts.kthermal        =  12;     % thermal conductivity
   opts.z_0             =  0.001;  % Surface aero. roughness length    [m]
   opts.z_obs           =  3.0;    % Height of wind and air temp obs   [m]
   opts.kabs_user       =  true;
   opts.tlagcolumn      =  6*3600/opts.dt;
   opts.use_ro_glc      =  false;
   opts.use_init        =  false;
   
   % define the thermal and spectral grid
   opts.dz_thermal      =  0.04;       % dz for heat transfer          [m]
   opts.dz_spectral     =  0.002;      % dz for radiative transfer     [m]
   opts.z0_thermal      =  20;         % thickness for heat transfer   [m]
   opts.z0_spectral     =  4;          % thickness for rad transfer    [m]
   opts.f_ice_min       =  0.01;
   
   opts.ro_snow_i       =  900.0;
   opts.i_grainradius   =  25;         % index 25 = 2.0 mm             [#]
   opts.nwavl           =  118;
   opts.nradii          =  35;
   
   % solver options
   opts.fzero           =  optimset('Display','off','TolX',1e-6);
   
   % model structure
   if strcmp(meltmodel,'icemodel')
      opts.skinmodel = false;
   elseif strcmp(meltmodel,'skinmodel')
      opts.skinmodel = true;
   end
   
   opts.simyears  = startyear:1:endyear;
   opts.numyears  = endyear-startyear+1;
