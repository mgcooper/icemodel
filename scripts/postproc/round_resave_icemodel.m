clean

% at the end of the day, ice2 dominates the storage, and it is already
% rounded. I might be able to save all mass fluxes in mm as integers

% I should also probably have not run the Leverett catchment and instead
% run 'region' for 2009-2011, then extracted Leverett. I cannot think of
% any reason not to do this. i could then also run with MODIS albedo
% instead of MAR, with SkinModel also. 

% note - for the skinmodel runs, i rounded the enbal vars from
% 1:length(fields)-3, meaning I rounded "Tflag". shouldn't matter, but I
% adjusted to 1:length(fields)-4 from thereafter

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% skinmodel
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% p.data  = ['/Volumes/Samsung_T5/matlab/GREENLAND/runoff/icemodel/model/' ...
%             'experiment/v8/region/skinmodel/output/reference/'];
% years   = 2009:2018;
% npts    = 551;
% 
% for k = 1:length(years)
%     year_k      = num2str(years(k));
%     path_k      = [p.data year_k '/'];
%     
%     for j = 1:npts
%         
%         % enbal
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
%         f       =   [path_k 'enbal/enbal_' int2str(j) '.mat'];
%         load(f)
%         fields   =   fieldnames(enbal);
%         
%         for i = 1:length(fields)-3
%             if strcmp(fields,'chi')
%                 enbal.(fields{i})   = roundn(enbal.(fields{i}),-3);
%             else        
%                 enbal.(fields{i})   = roundn(enbal.(fields{i}),-1);
%             end
%         end
%         save(f,'enbal');
%         
%         f       =   [path_k 'ice1/ice1_' int2str(j) '.mat'];
%         load(f)
%         
%         % ice1
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
%         fields   =   fieldnames(ice1);
%         for i = 1:length(fields)-3
%             if strcmp(fields,'Tice')
%                 ice1.(fields{i}) = roundn(ice1.(fields{i}),-1);
%             else        
%                 ice1.(fields{i}) = roundn(ice1.(fields{i}),-4);
%             end
%         end
%         save(f,'ice1');
%     end
% end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% icemodel
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p.data  = ['/Volumes/Samsung_T5/matlab/GREENLAND/runoff/icemodel/model/' ...
            'experiment/v8/region/icemodel/output/'];
years   = 2012:2018;
npts    = 1487;
folds   = {'reference/','t1/','t2/','t3/','t4/'};
% m=3, k =3, j=164, so run j=164 to end, then k = 4:end, then restart m=4
for m = 4:length(folds)
    path_m          = [p.data folds{m}];

    for k = 1:length(years)
        year_k      = num2str(years(k));
        path_k      = [path_m year_k '/'];

        for j = 1:npts

            % start enbal
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            f       =   [path_k 'enbal/enbal_' int2str(j) '.mat'];
            load(f)
            fields   =   fieldnames(enbal);

            for i = 1:length(fields)-4
                if strcmp(fields,'chi')
                    enbal.(fields{i})   = roundn(enbal.(fields{i}),-3);
                else        
                    enbal.(fields{i})   = roundn(enbal.(fields{i}),-1);
                end
            end
            save(f,'enbal');
            % end enbal
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            % start ice 1
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            f       =   [path_k 'ice1/ice1_' int2str(j) '.mat'];
            load(f)

            fields   =   fieldnames(ice1);
            for i = 1:length(fields)-3
                if strcmp(fields,'Tice')
                    ice1.(fields{i}) = roundn(ice1.(fields{i}),-1);
                else        
                    ice1.(fields{i}) = roundn(ice1.(fields{i}),-4);
                end
            end
            save(f,'ice1');
            % end ice 1
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
        end
    end
end

        