function [Invertf,sigmaP,bP,RelTotalError,AbsTotalError,AbsTruncateError,AbsRoundoffError,LaguerreCoef] =...
wfnWeeksCoreSigmab(varargin)
%WFNWEEKSCORESIGMAB The core function for the Weeks' method computations
%  for a single time with \sigma and b as parameters.
%  The outer script selects the appropriate parameters and sends them to this core function. 
%
%  Use:
%  [Invertf,sigmaP,bP,RelTotalError,AbsTotalError,AbsTruncateError,AbsRoundoffError,LaguerreCoef] =...
%  wfnWeeksCoreSigmab(varargin)
%  
%  Input:
%  (variable inputs)
%  1) FLaplace = a symbolic expression for the Laplace transform space function F(s)
%  2) TimeInput = a single time where f(t) is to be computed
%  3) NLag = the number of terms in the Laguerre expansion
%  if 5 inputs
%   4) sigmaP = user defined free parameter sigma in the Weeks method  
%   5) bP = user defined free parameter b in the Weeks method 
%  if 9 inputs
%   4) sigmin = the minimum sigma value for the (sigma,b) search
%   5) sigmax = the maximum sigma value for the (sigma,b) search
%   6) bmax = the maximum b value for the (sigma,b) search
%   7) SearchSwitch = 0-CPU local fminbnd, 1-CPU global
%   8) tols = the resolution of the s-parameter
%   9) tolb = the resolution of the b-parameter 
%
%  Output: 
%  Invertf = numerically computed f(t)
%  sigmaP = automatically determined sigma or user defined if an input
%  bP = automatically determined b or user defined if an input
%  RelTotalError = the estimate relative total error
%  AbsTotalError = the estimated absolute total error 
%  AbsTruncationError = the estimated absolute trunction error 
%  AbsRoundoffError = the estimated absolute round-off error 
%  LaguerreCoef = Laguerre polynomial expansion coefficients
%
%  Author: 
%  Patrick Kano, Moysey Brio - 2011, 2016
%
%  Modification Date [M/D/Y]:
%  03/31/2011 - Initial work
%  04/05/2016 - Jacket removed and replaced with a (rho,alpha) 2D search

switch nargin
 case 5 %input sigma and b parameters input by the user
  FLaplace = varargin{1};
  TimeInput = varargin{2};
  NLag = varargin{3};
  sigmaP = varargin{4};
  bP = varargin{5};
 
 case 9 %determine the optimal sigma and b parmeters internally
  FLaplace = varargin{1};
  TimeInput = varargin{2};
  NLag = varargin{3};
  sigmin = varargin{4};
  sigmax = varargin{5};
  bmax = varargin{6}; 
  SearchSwitch = varargin{7};
  tols = varargin{8};
  tolb = varargin{9};
  
  [sigmaP, bP] = wfnParamEstSigmab(FLaplace,TimeInput,NLag,sigmin,sigmax,bmax,SearchSwitch,tols,tolb);
  
 otherwise
  error('wfnWeeksCoreSigmab','Incorrect number of inputs, 5 or 9 are required.');
end
 
%Laguerre polynomial coefficients calculation 
LaguerreCoef = wfncpuFFTLagCoefSigmab(FLaplace,NLag,sigmaP,bP);

%Clenshaw algorithm to sum the polynomials
Invertf = wfnClenshawSigmab(NLag,TimeInput,sigmaP,bP,LaguerreCoef);
 
%Estimate the absolute and relative error of the inverted function 
%Note: Reversing order of b,sigma from previous two functions is correct.
[AbsTotalError, AbsTruncateError, AbsRoundoffError] = wfncpuErrorEstSigmab(bP,sigmaP,FLaplace,TimeInput,NLag,0); 

%Convert from the log base 10 value
AbsTotalError = 10^(AbsTotalError);
AbsTruncateError = 10^(AbsTruncateError);
AbsRoundoffError = 10^(AbsRoundoffError);

RelTotalError = abs(AbsTotalError)/abs(Invertf);

end %function definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
