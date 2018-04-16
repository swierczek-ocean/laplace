function TotalError = wfncpuErrorEstVecAR(yvec,FLaplace,TimeInput,NLag,ErrorFlag)
%WFNCPUERRORESTVECAR A function to return the absolute error estimate of the approximate f(t) from the Weeks method
%  The function returns the total, truncation, and round-off error 
%  from the tail of the sequence of Laguerre coefficients.
%  This function uses standard MATLAB/CPU computations. 
%
%  Use: 
%  TotalError = wfncpuErrorEstVecAR(yvec,FLaplace,TimeInput,NLag,ErrorFlag)
%  
%  Input:
%  yvec = [rhoP,alphaP]
%  rhoP = Weeks rho parameter
%  alphaP = Weeks alpha parameter
%  FLaplace = a symbolic expression for the Laplace transform space function F(s)
%  TimeInput = time where f(t) is to be computed
%  NLag = number of terms in the Weeks expansion
%  ErrorFlag = flag to define the error type
%   0 = log base 10 of the total error
%   1 = log base 10 of the truncation error only
%
%  Output:
%  TotalError = total/sum of absolute truncation and round-off error estimates
%  TruncateError = absolute truncation error estimate
%  RoundoffError = absolute roundoff error estimate 
%
%  Author: 
%  Patrick Kano, Moysey Brio - 2016
%
%  Modification Date [M/D/Y]:
%  04/09/2016 - Adapationl of the error estimate fuction to take a vector
%  input for alpha and rho

rhoP = yvec(1);
alphaP = yvec(2);

[TotalError, TruncateError, RoundoffError] =...
wfncpuErrorEstAlphaRho(rhoP,alphaP,FLaplace,TimeInput,NLag,ErrorFlag);

end %function definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
