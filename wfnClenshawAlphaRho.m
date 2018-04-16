function ftimevec = wfnClenshawAlphaRho(NLag,TimeInput,alphaP,rhoP,LaguerreCoef)
%WFNCLENSHAWALPHARHO Laguerre polynomial evaluation 
%  The function evaluates the Laguerre polynomial
%  expansion with coefficients {a} at 
%  the values {t} using Clenshaw's algorithm for the Weeks method.   
%  
%  Use:
%  ftimevec = wfnClenshawAlphaRho(NLag,TimeInput,alphaP,rhoP,LaguerreCoef)
%  
%  Input:
%  NLag = number of Laguerre polynomial coefficients
%  TimeInput = input time (scalar or a column vector)
%  alphaP = Weeks' method parameter alpha
%  rhoP = Week's method parameter rho
%  LaguerreCoef = Laguerre polynomial expansion coefficients
%
%  Output: 
%  ftimevec = the summed Laguerre polynomial expansion
%  
%  Comment:
%  The Clenshaw algorithm allows for a stable recursive evalation of
%  the Laguerre polynomials.   
%
%  References: 
%  Clenshaw algorithm
%  www.math.arizona.edu/~swig/documentation/matlab-tips/
%  www.en.wikipedia.org/wiki/Clenshaw_algorithm
%
%  Author: 
%  Patrick Kano, Moysey Brio - 2016
%
%  Modification Date [M/D/Y]:
%  03/04/2016 - Version 1.0

%Times are all independant
Ntimes = length(TimeInput);

Cpresent = zeros(Ntimes,1);
Cptwo = zeros(Ntimes,1);
Cpone = LaguerreCoef(NLag)*ones(Ntimes,1);

for kidx=(NLag-1):-1:1
 Cpresent = ((2*kidx-1-(rhoP*TimeInput))/kidx).*Cpone - (kidx/(kidx+1)).*Cptwo + LaguerreCoef(kidx)*ones(Ntimes,1);
 Cptwo = Cpone; 
 Cpone = Cpresent;  
end
 
ftimevec = exp((-1.0*alphaP)*TimeInput).*Cpresent;

end %function definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
