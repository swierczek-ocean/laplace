function TotalError = wfnNestedErrorSigmab(sigmaP,FLaplace,TimeInput,NLag,bmax,tolb)
%WFNNESTEDERROR The cost function in the local fminbnd serach for the (sigma,b) parameters.
%  For a fixed sigma, the truncation error is minimized with respect to the b parameter.
%  The total error is computed for the (sigma,b) pair and minimized.
% 
%  Use: 
%  [TotalError, TruncateError, RoundoffError] = wfnNestedError(sigmaP,FLaplace,TimeInput,NLag,bmax,tolb)
%  
%  Input:
%  sigmaP = Weeks sigma parameter
%  FLaplace = a symbolic expression for the Laplace transform space function F(s)
%  TimeInput = ordinate where f(t) are to be computed (scalar).
%  NLag = number of terms in the Weeks expansion.
%  bmax = maximum possible value of b.
%  tolb = 'TolX' option to fminbnd for the mininimzation of the truncation error over bP
%
%  Output:
%  TotalError = total absolute error estimate 
%  
%  Author: 
%  Patrick Kano, Moysey Brio - 2011
%
%  Modification Date [M/D/Y]:
%  04/01/2011 - Version 1.0

ErrorFlag = 1; %Truncation error only
options = optimset('TolX',tolb);
bP = fminbnd('wfncpuErrorEstSigmab',0,bmax,options,sigmaP,FLaplace,TimeInput,NLag,ErrorFlag);

ErrorFlag = 0; %Total error
TotalError = wfncpuErrorEstSigmab(bP,sigmaP,FLaplace,TimeInput,NLag,ErrorFlag);

end %function definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
