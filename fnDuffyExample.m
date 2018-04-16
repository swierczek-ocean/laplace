function Fs = fnDuffyExample(s,rho,r)
%FNDUFFYEXAMPLE Duffy Numerical Methods F(s) Example
%   An implementation of the Laplace space function F(s) in the
%   book by Dean Duffy 
%   Green's Functions with Applications
%   Equation 7.2.20
%  
%   April 8, 2016
%
%   G(r,s|rho,0^{+}) = \frac{\pi
%   \rho}{s^{1/3}}*\{H(r-\rho)*Bi(\eta_{*})*Ai(\eta) +
%   [1-H(r-\rho)]*Ai(\eta_{*})*Bi(\eta) - C*Ai(\eta_{*})*Ai(\eta)

%Airy functions Ai and Modified Bessel Function Bi

%Parameters
%rho = 3;
%r = 3;

%Laplace space function
eta     = (s.^(1/3)).*(r + 1./(4*s));
etastar = (s.^(1/3)).*(rho + 1./(4*s));

Besseletastar = besselj(1,etastar);
Besseleta     = besselj(1,eta);
Airyetastar   = airy(1,etastar);
Airyeta       = airy(eta);

Kappa = (pi*rho*exp(0.5*(rho-r)))./(s.^(1/3));

Ccoeff = 1; %What should this be?

if(r==rho)
 Fs = 0.5*Kappa.*(Besseletastar.*Airyeta - Ccoeff.*Airyetastar.*Airyeta)+...
 0.5*Kappa.*(Besseleta.*Airyetastar - Ccoeff.*Airyetastar.*Airyeta);   
elseif(rho<r) %r-rho > 0
 Fs = Kappa.*(Besseletastar.*Airyeta - Ccoeff.*Airyetastar.*Airyeta);    
else
 Fs = Kappa.*(Besseleta.*Airyetastar - Ccoeff.*Airyetastar.*Airyeta);   
end

end %function definitio

