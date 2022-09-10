%--------------------------------------------------------------------------
%   COMPUTE THE STABILITY EXCHANGE COEFFICIENTS
%--------------------------------------------------------------------------

function [D_h,D_e] = EXCOEFS_SCALAR(z_0,z_windobs,windspd,xkappa)

% this is matt hoffman's version. it isn't implemented yet.

	visc_air = 1.461e-5;
	alpha_h =  1.0;
	alpha_q =  1.0;

% Note: u_fric should account for stability!!!!!
	u_fric = xkappa*z_0 / (log(z_windobs/z_0));
	Re = u_fric*z_0/visc_air;
	
    if (Re<0.135)
		z_t = z_0*exp(1.250);
		z_q = z_0*exp(1.610);
	elseif (Re<2.5)
		z_t = z_0*exp(0.149-0.55*log(Re));
		z_q = z_0*exp(0.351-0.628*log(Re));
	else
		z_t = z_0*exp(0.317-0.565*log(Re)-0.183*(log(Re))^2);
		z_q = z_0*exp(0.396-0.512*log(Re)-0.180*(log(Re))^2);
    end
    
	CD = xkappa^2/(log(z_windobs/z_0))^2;
	CH = alpha_h*xkappa*sqrt(CD)/(xkappa*CD^(-0.5) - log(z_t/z_0));
	CE = alpha_q*xkappa*sqrt(CD)/(xkappa*CD^(-0.5) - log(z_q/z_0));
	D_h = CH * windspd;
	D_e = CE * windspd;
    
%	Qh = ro_air*Cp*CH*wspd*stability*(Tair-T_sfc)
%	Qe = ro_air*xLs*CE*wspd*stability*CCC*(ea-es0)
    