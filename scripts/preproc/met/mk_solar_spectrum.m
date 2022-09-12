clean

fname               =   'solar.dat'; % output file name
save_data           =   1;
plot_data           =   0;
%%

homepath            =   pwd;

if strcmp(homepath(2:6),'Users')
    path.data       =   ['/Users/mattcooper/Dropbox/CODE/MATLAB/GREENLAND/' ...
                            'icemodel/data/input/'];
    path.save       =   ['/Users/mattcooper/Dropbox/CODE/MATLAB/GREENLAND/' ...
                            'icemodel/data/output/'];
elseif strcmp(homepath(10:16),'mcooper')    
end

%%

% read in my solar spectrum
load([path.data 'Qi_coart.mat']);
QiCoart             =   Qi_coart.Qi;
lambdaCoart         =   Qi_coart.Lambda;

% read in Glen's spectrum
solar               =   importdata([path.data 'solar.dat']);
lambda              =   1000.*solar(:,1);
QiListon            =   solar(:,2);

% interpolate my spectrum to his lambda
Qi                  =   interp1(lambdaCoart,QiCoart,lambda);

% fill in missing lambda with his values 
Qi(isnan(Qi))       =   QiListon(isnan(Qi));


%% make figure

if plot_data == 1
    figure(1)
    p(1)            =   plot(lambda,solar(:,2)); hold on;
    p(2)            =   plot(lambdaCoart,QiCoart);
    p(3)            =   plot(lambda,Qi);
    l               =   legend('Antarctica','Greenland - OG','Greenland - Interp');
    if save_fig == 1
        export_fig([path.save 'Qi_incoming.png'],'-r400')
    end
end

%% prep data for output

data                =   roundn([lambda./1000 Qi],-5);
data                =   data';
pLambda             =   '%8.5f';
pQi                 =   '%12.5f\n';

%% save the output

if save_data == 1
    
    fid             =   fopen([path.save fname],'w'); 
    fprintf(fid,[pLambda pQi],data);
    fclose(fid);
end