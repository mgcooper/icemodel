function spect_coefs = SPECTEXTCOEF(opts,qext,g,ss_coalb,r_snow)
%SPECTEXTCOEF compute the spectral extinction coefficients

%  bulk ext coeffs are updated each timestep with changing ice density
   iradius        =  opts.i_grainradius;
   sigma_e        =  0.75.*qext(iradius,:)./r_snow;
   spect_coefs    =  sigma_e .* sqrt(ss_coalb(iradius,:) - ...
                     ss_coalb(iradius,:) .* g(iradius,:) + ...
                     ss_coalb(iradius,:).^2 .* g(iradius,:));
