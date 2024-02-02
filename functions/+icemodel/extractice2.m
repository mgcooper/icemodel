function vmat = extractice2(ice2, varname, t1, t2, z1, z2, dz, dt)
   %EXTRACTICE2 Extract 2-d icemodel data

   arguments
      ice2
      varname % = inputname(1)
      t1 = 1
      t2 = nan
      z1 = 0
      z2 = 2
      dz = 0.04
      dt = 3600; % 900 use 3600 b/c this more often called after POSTPROC
   end

   if isdatetime(t1) && isdatetime(t2)

   end


   tidx = t1:t2;
   zidx = max(1, floor(z1 / dz(1))) : ceil(z2 / dz(1));
   vmat = ice2.(varname)(zidx, tidx);
   tmat = dt * (t1:t2);
   zmat = 0:dz(1):z2;

   if numel(zmat) ~= size(vmat, 1)
      zmat = 0:dz(1):z2-dz(1);
   end


   % input processing from plotice2, not developed for this function
   % if isnumeric(varname)
   %    varname = inputname(1);
   % end
   %
   % if nargin < 4
   %    z2 = 2;
   % end
   % if nargin < 5
   %    try
   %       dz = ice2.Z(2) - ice2.Z(1);
   %    catch
   %       dz = 0.04;
   %    end
   % end
   % if nargin < 6
   %    dt = 3600;
   % end

end
