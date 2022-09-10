function out = linearfill(in)

% in = the data you want to fill
% si = starting index
% ini = starting value (value of in at index si)
% 
for n = 1:length(in)
    
    if isnan(in(n))
        
        si              =   n;
        if n > 1
            ini         =   in(n-1);
        else
            ini         =   in(n);
        end
        
        i               =   1;
        ii              =   1;
        
        while ii == i
            if n+i > length(in)             % if we reach the end of in
                ii      =   ii-2;
                ei      =   n+ii+1;
                in(ei)  =   ini;
                break
            end
            if isnan(in(n+i))
                i       =   i+1;
                ii      =   ii+1;
            else
                ii      =   ii-1;
            end
        end
        
        ei              =   n+ii+1;
        ine             =   in(ei);
        if n > 1
            infill      =   cumsum(((ine-ini)/(ei-si)).*ones((ei-si),1));
            infill      =   in(si-1)+infill;
        else
            infill      =   ine;
        end
        
        in(si:ei-1)     =   infill;
        
    end
end

out = in;

