function uhn = wfnTRACIcoreUH(theta,Fofs,alphaP,rhoP,Nmin,Nmax)
%WFNTRACICOREUH Core TRACI function for UH_{n}
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

Hoftheta = real(fcore);

wterm = cos(theta*(Nmin:Nmax))';

uhn = real((1/(2*pi))*( wterm.*Hoftheta ));

end %function definition

