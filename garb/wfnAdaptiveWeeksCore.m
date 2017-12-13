function aintegrand = wfnAdaptiveWeeksCore(w,Fofs,alphaP,rhoP,Nmin,Nmax)
%WFNADAPTIVEWEEKSCORE Core complex function integrated in the Weeks method
%   This function is called by the adaptive integration fuction when
%   computing the coefficients in the Weeks method
%   Note that this is a vectorized implemetation that returns the
%   core functios from Nmin to Nmax.
%   Patrick Kano, Moysey Brio
%   June 6, 2016

s = rhoP./(1-w) - alphaP;

Fofw = Fofs(s);

fcore = (1/(2i*pi))*((rhoP./(1-w)).*Fofw);

if isscalar(s)
 aintegrand = (w.^(-1.0*(Nmin:Nmax))).*fcore;
else
 wterm = w.^(-1.0*(Nmin:Nmax)');
 aintegrand =  wterm.*fcore;  
     
 %for jidx=1:length(s)
 % aintegrand(:,jidx) = (w.^(-1.0*(Nmin:Nmax))).*fcore(jidx);
 %end %jidx
end %if-else
end %function definition

