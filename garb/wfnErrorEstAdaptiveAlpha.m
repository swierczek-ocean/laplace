function [TotalError, TruncateError, RoundoffError] = wfnErrorEstAdaptiveAlpha(alphaP,FLaplace,NLag,ErrorFlag)
%WFNCPUERRORESTADAPTIVEALPHA A function to return the absolute error estimate of the approximate f(t) from the Weeks method
%  when one uses adaptive integration to estimate the Laguerre expansion coefficients.
%  The function returns the total, truncation, and round-off error 
%  from the tail of the sequence of Laguerre coefficients.
%  This function uses the MATLAB adaptive integration approach with waypoints
%  along a unit circle in the W- complex plane.
%
%  Use: 
%  [TotalError, TruncateError, RoundoffError] = wfnErrorEstAdaptiveAlpha(alphaP,FLaplace,NLag,ErrorFlag)
%  
%  Input:
%  alphaP = Weeks alpha parameter
%  FLaplace = a symbolic expression for the Laplace transform space function F(s)
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
%  06/10/2016 - Initial Implementation

%Double the number of coefficients actually used to get the tail of the sequence
Ntail = 2*NLag; 

rhoP = 1.0;

%RelTol and AbsTol could be adjusted
%Waypoints could also be adjusted

if(1)
 iNmin = 1;
 iNmax = Ntail;
 
 av =...
 integral(@(w)wfnAdaptiveWeeksCore(w,FLaplace,alphaP,rhoP,iNmin,iNmax),1,1,...
 'ArrayValued',true,...
 'Waypoints',[1i,-1,-1i],...
 'RelTol',0.001,...
 'AbsTol',0.001);
else
 tmpWayPoints = [(sqrt(2)/2)*(1+1i),1i,(sqrt(2)/2)*(-1+1i),-1,(sqrt(2)/2)*(-1-1i),-1i,(sqrt(2)/2)*(1-1i)];
 for nidx=1:Ntail    
  av(nidx) =....
  quadgk(@(w)wfnAdaptiveWeeksCore(w,FLaplace,alphaP,rhoP,nidx,nidx),1,1,...
 'Waypoints',tmpWayPoints,...
 'AbsTol',0.01,...
 'RelTol',0.01);
 end %nidx
end

% ErrorFlag==0 -> Total error, ErrorFlag==1 -> Truncation error only
RoundoffError = (1-ErrorFlag)*eps*sum( abs( av(1:NLag) ) );
TruncateError = sum( abs( av(NLag+1:2*NLag) ) );

%I do not include time in the estimate
%TotalError = (1-ErrorFlag)*(exp(-alphaP*TimeInput)*(TruncateError+RoundoffError)) + ErrorFlag*TruncateError;
TotalError = (1-ErrorFlag)*((TruncateError+RoundoffError)) + ErrorFlag*TruncateError;

RoundoffError = log10(RoundoffError);
TruncateError = log10(TruncateError);
TotalError = log10(TotalError);

end %function definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
