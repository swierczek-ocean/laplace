function [TotalError, TruncateError, RoundoffError] = wfncpuErrorEstSigmab(bP,sigmaP,FLaplace,TimeInput,NLag,ErrorFlag)
%WFNCPUERRORESTSIGMAB A function to return the absolute error estimate of the approximate f(t) from the Weeks method
%  The function returns the total, truncation, and round-off error 
%  from the tail of the sequence of Laguerre coefficients.
%  This function uses standard MATLAB/CPU computations. 
%
%  Use: 
%  [TotalError, TruncateError, RoundoffError] = wfncpuErrorEstSigmab(bP,sigmaP,FLaplace,TimeInput,NLag,ErrorFlag)
%  
%  Input:
%  bP = Weeks b parameter
%  sigmaP = Weeks sigma parameter
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
%  Patrick Kano, Moysey Brio - 2011
%
%  Modification Date [M/D/Y]:
%  03/20/2011 - Version 1.0
%  03/04/2016 - Simply a name change to signify the (sigma, b) parameters

%The only only difference betwen cpu and jacket versions of this function is how the FFT is performed.
%Separate CPU/GPU FFT functions used to avoid logical branches.

Ntail = 2*NLag; %Double the number of coefficients actually used to get the tail of the sequence

av = wfncpuFFTLagCoefSigmab(FLaplace,Ntail,sigmaP,bP);

% ErrorFlag==0 -> Total error, ErrorFlag==1 -> Truncation error only
RoundoffError = (1-ErrorFlag)*eps*sum( abs( av(1:NLag) ) );
TruncateError = sum( abs( av(NLag+1:2*NLag) ) );
TotalError = (1-ErrorFlag)*(exp(sigmaP*TimeInput)*(TruncateError+RoundoffError)) + ErrorFlag*TruncateError;

RoundoffError = log10(RoundoffError);
TruncateError = log10(TruncateError);
TotalError = log10(TotalError);

end %function definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
