function [Invertf,alphaP,rhoP,RelTotalError,AbsTotalError,AbsTruncateError,AbsRoundoffError,LaguerreCoef] =...
wfnWeeksCoreAdaptiveIntegrate(varargin)
%WFNWEEKSCOREADAPTIVEINTEGRATE The core function for the Weeks' method computations.
%  Two modalities are possible depending on the number of input parameters.
%  Modality 1: for a single alpha, return f(t)
%  Modality 2: for a range of alpha values, estimate the alpha value that
%  yields the best approximation for f(t).  The use it to estimate f(t)
%
%  Use:
%  [Invertf,alphaP,rhoP,RelTotalError,AbsTotalError,AbsTruncateError,AbsRoundoffError,LaguerreCoef] =...
%  wfnWeeksCoreAdaptiveIntegrate(varargin)
%  
%  Input:
%  (variable inputs)
%  1) FLaplace = a symbolic expression for the Laplace transform space function F(s)
%  2) TimeInput = a single time where f(t) is to be computed
%  3) NLag = the number of terms in the Laguerre expansion
%  if 4 inputs
%   4) alphaP = user defined free parameter alpha in the Weeks method  
%  if 6 inputs
%   4) alphamin = the minimum alpha value for the search
%   5) alphamax = the maximum alpha value for the search
%   6) tolalpha = the resolution of the alpha-parameter
%
%  Output: 
%  Invertf = numerically computed f(t)
%  alphaP = automatically determined alpha or user defined if an input
%  rhoP = automatically determined rho or user defined if an input
%  RelTotalError = the estimate relative total error
%  AbsTotalError = the estimated absolute total error 
%  AbsTruncationError = the estimated absolute trunction error 
%  AbsRoundoffError = the estimated absolute round-off error 
%  LaguerreCoef = Laguerre polynomial expansion coefficients
%
%  Author: 
%  Patrick Kano, Moysey Brio - 2016
%
%  Modification Date [M/D/Y]:
%  06/10/2016 - Initial version with alpha search and rho = 1

switch nargin
 case 4 %input alpha parameters from the user
  FLaplace = varargin{1};
  TimeInput = varargin{2};
  NLag = varargin{3};
  alphaP = varargin{4};

  rhoP = 1.0;
 case 6 %determine the optimal alpha
  FLaplace = varargin{1};
  TimeInput = varargin{2}; %Not used for the estimated optimal parameters
  NLag = varargin{3};
  alphamin = varargin{4};
  alphamax = varargin{5};
  tolalpha = varargin{6};

  [alphaP, rhoP] =...
  wfnParamEstAdaptiveAlpha(FLaplace,NLag,alphamin,alphamax,tolalpha);
  
 otherwise
  error([mfilename,':nargin'],'Incorrect number of inputs, 4 or 6 are required.');
end

%Laguerre polynomial coefficients calculation
if(1)
 iNmin = 1; %why not 0 to NLAG?
 iNmax = NLag;
 
 LaguerreCoef =...
 integral(@(w)wfnAdaptiveWeeksCore(w,FLaplace,alphaP,rhoP,iNmin,iNmax),1,1,...
 'ArrayValued',true,...
 'Waypoints',[1i,-1,-1i],...
 'RelTol',0.01,...
 'AbsTol',0.01);
else
 %tmpWayPoints = [(sqrt(2)/2)*(1+1i),1i,(sqrt(2)/2)*(-1+1i),-1,(sqrt(2)/2)*(-1-1i),-1i,(sqrt(2)/2)*(1-1i)];
 for nidx=1:NLag
  LaguerreCoef(nidx) =....
  quadgk(@(w)wfnAdaptiveWeeksCore(w,FLaplace,alphaP,rhoP,nidx,nidx),1,1,...
 'Waypoints',[1i,-1,-1i],...
 'AbsTol',0.01,...
 'RelTol',0.01);
 end %nidx
end
LaguerreCoef
if(norm(LaguerreCoef)<1e-18) %quadrature failure
 Invertf = NaN
 
 AbsTotalError = Inf;
 AbsTruncateError =  Inf;
 AbsRoundoffError =  Inf;

 RelTotalError =  Inf;
 
else
 %Clenshaw algorithm to sum the polynomials
 Invertf = wfnClenshawAlphaRho(NLag,TimeInput,alphaP,rhoP,LaguerreCoef);
 
 %Estimate the absolute and relative error of the inverted function 
 %Note: Reversing order of rho,alpha for this function is correct.
 ErrorFlag = 0;
 [AbsTotalError, AbsTruncateError, AbsRoundoffError] =...
 wfnErrorEstAdaptiveAlpha(alphaP,FLaplace,NLag,ErrorFlag); 

 %Convert from the log base 10 value
 AbsTotalError = 10^(AbsTotalError);
 AbsTruncateError = 10^(AbsTruncateError);
 AbsRoundoffError = 10^(AbsRoundoffError);

 RelTotalError = abs(AbsTotalError)/abs(Invertf);
end
%This plots the trasformations of the contour (optional)
%fnSequenceTransforms(alphaP,rhoP,101)

end %function definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
