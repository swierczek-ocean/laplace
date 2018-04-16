function TotalError = wfnNestedErrorAlphaRho(alphaP,FLaplace,TimeInput,NLag,rhomax,tolrho)
%WFNNESTEDERRORALPHARHO The cost function in the local fminbnd serach for the (alpha,rho) parameters.
%  For a fixed alpha, the truncation error is minimized with respect to the rho parameter.
%  The total error is computed for the (alpha,rho) pair and minimized.
% 
%  Use: 
%  [TotalError, TruncateError, RoundoffError] = wfnNestedError(alphaP,FLaplace,TimeInput,NLag,rhomax,tolrho)
%  
%  Input:
%  alphaP = Weeks alpha parameter
%  FLaplace = a symbolic expression for the Laplace transform space function F(s)
%  TimeInput = ordinate where f(t) are to be computed (scalar).
%  NLag = number of terms in the Weeks expansion.
%  rhomax = maximum possible value of rho.
%  tolrho = 'TolX' option to fminbnd for the mininimzation of the truncation error over rhoP
%
%  Output:
%  TotalError = total absolute error estimate 
%  
%  Author: 
%  Patrick Kano, Moysey Brio - 2016
%
%  Modification Date [M/D/Y]:
%  04/08/2016 - Version 1.0, Usig total error instead of just truncation error

ErrorFlag = 0; %Total error
options = optimset('TolX',tolrho);
rhoP = fminbnd('wfncpuErrorEstAlphaRho',0,rhomax,options,alphaP,FLaplace,TimeInput,NLag,ErrorFlag);

ErrorFlag = 0; %Total error
TotalError = wfncpuErrorEstAlphaRho(rhoP,alphaP,FLaplace,TimeInput,NLag,ErrorFlag);

end %function definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
