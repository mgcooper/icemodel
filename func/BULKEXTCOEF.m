%--------------------------------------------------------------------------
%   Compute the spectrally integrated extinction coefficient. 
%--------------------------------------------------------------------------
function bulkcoefs = BULKEXTCOEF(dz_spect,ro_sno,spect_lower,spect_upper,solardwavl)
%--------------------------------------------------------------------------

% Compute the downward bulk extinction coefficient, scaled by the total ice
%   equivalent thickness of each layer
   bulkcoefs = -log((sum(exp(spect_lower.*ro_sno).*solardwavl,2))./    ...
         (sum(exp(spect_upper.*fix(max(ro_sno,300))).*solardwavl,2)))./dz_spect;

% Cast in a form that can be used in the two-stream computation (add the
%   boundaries).  Here I have assumed that it is okay to call the
%   boundaries equal to the value at the center of that grid cell (the
%   value prevails thoughout the cell).
   bulkcoefs = vertcat(bulkcoefs,bulkcoefs(end),bulkcoefs(end));
   
   
% % to run these tests, uncomment but first run the 

% % keep a copy before multiplying by ro_sno/ro_ice
%    speckeep    = speccoefs;

% % Explicit, but slow version:
%    ro_sno      = max(ro_sno,300);
%    zj          = z_walls(1:end-1);
%    zjj         = z_walls(2:end);
%    dz          = zjj(2)-zjj(1);
%    speccoefs   = (ro_sno./917).*speccoefs;
%    sum1        = sum(exp(-speccoefs.*zjj).*solardwavl);
%    sum2        = sum(exp(-speccoefs.*zj).*solardwavl);
%    bulkcoefs   = -log(sum1./sum2)./dz;
   
% % here i test using a pre-computed solar.*exp(-spec_coefs.*z)
%    spec_walls  = exp(-z_walls.*repmat(speckeep./917,1,2001));
%    spec_walls1 = spec_walls(:,2:end);
%    spec_walls2 = spec_walls(:,1:end-1);
%    solardwavl  = solar.*dwavl;
%    test1       = sum(solardwavl.*spec_walls1.^ro_sno);
%    test2       = sum(solardwavl.*spec_walls2.^ro_sno);
%    bulk_test   = -log(test1./test2)./dz;
%    
%    figure; myscatter(test1,sum1); addOnetoOne;
%    figure; myscatter(test2,sum2); addOnetoOne;
% 
% % this tests not taking exp() first
%    spec_walls  = -z_walls.*repmat(speckeep./917,1,2001);
%    spec_walls1 = spec_walls(:,2:end);
%    spec_walls2 = spec_walls(:,1:end-1);
%    test1       = sum(exp(spec_walls1.*ro_sno).*solardwavl);
%    test2       = sum(exp(spec_walls2.*ro_sno).*solardwavl);
%    bulk_test   = -log(test1./test2)./dz;
%    
%    figure; myscatter(test1,sum1); addOnetoOne;
%    figure; myscatter(test2,sum2); addOnetoOne;
   
% % here i was testing using the integral solution (evaluate endponts)
%    test1       = test(:,2:end).^(ro_snow./ro_ice).*(solar(end)-solar(1));
%    test2       = solar.*test(:,1:end-1).^(ro_snow./ro_ice).*dwavl);
%    
%    sum1        = -solar./spec_coefs.*exp(-spec_coefs.*zjj);
%    sum2        = -solar./spec_coefs.*exp(-spec_coefs.*zj);

   
