function ujn = wfnTRACIcoreUJ(theta,Fofs,alphaP,rhoP,Nmin,Nmax)
%WFNTRACICOREUJ Core TRACI function for UJ_{n}
%   This function is called by the adaptive integration fuction when
%   computing the coefficients in the Weeks method
%   Note that this is a vectorized implemetation that returns the
%   core functions from Nmin to Nmax.
%   Patrick Kano, Moysey Brio
%   June 22, 2016

%   TRACI - Weeks Method Laplace Transform by Adaptively Computed Integrals
%   fcore = H + iJ
cost = cos(theta);
sint = sin(theta);
s = rhoP./(1-cost-1i*sint) - alphaP;

fcore = (rhoP./(1-cost-1i*sint)).*Fofs(s);

Joftheta = imag(fcore);

wterm = sin(theta*(Nmin:Nmax))';

ujn = real((1/(2*pi))*( wterm.*Joftheta ));

end %function definition

