function [Invertf,alphaP,rhoP,RelTotalError,AbsTotalError,AbsTruncateError,AbsRoundoffError,LaguerreCoef] =...
wfnWeeksCoreAlphaRho(varargin)
%WFNWEEKSCOREALPHARHO The core function for the Weeks' method computations
%  for a single time with \alpha and \rho as parameters.
%  The outer script selects the appropriate parameters and sends them to this core function. 
%
%  Use:
%  [Invertf,alphaP,rhoP,RelTotalError,AbsTotalError,AbsTruncateError,AbsRoundoffError,LaguerreCoef] =...
%  wfnWeeksCoreAlphaRho(varargin)
%  
%  Input:
%  (variable inputs)
%  1) FLaplace = a symbolic expression for the Laplace transform space function F(s)
%  2) TimeInput = a single time where f(t) is to be computed
%  3) NLag = the number of terms in the Laguerre expansion
%  if 5 inputs
%   4) alphaP = user defined free parameter alpha in the Weeks method  
%   5) rhoP = user defined free parameter rho in the Weeks method 
%  if 9 inputs
%   4) alphamin = the minimum alpha value for the (alpha,rho) search
%   5) alphamax = the maximum alpha value for the (alpha,rho) search
%   6) rhomax = the maximum rho value for the (alpha,rho) search
%   7) SearchSwitch = search switch passed to wfnParamEstAlphaRho
%   8) tolalpha = the resolution of the alpha-parameter
%   9) tolrho = the resolution of the rho-parameter 
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
%  03/04/2016 - Initial version with alpha,rho searching

switch nargin
 case 5 %input sigma and b parameters input by the user
  FLaplace = varargin{1};
  TimeInput = varargin{2};
  NLag = varargin{3};
  alphaP = varargin{4};
  rhoP = varargin{5};
 
 case 9 %determine the optimal sigma and b parmeters internally
  FLaplace = varargin{1};
  TimeInput = varargin{2};
  NLag = varargin{3};
  alphamin = varargin{4};
  alphamax = varargin{5};
  rhomax = varargin{6}; 
  SearchSwitch = varargin{7};
  tolalpha = varargin{8};
  tolrho = varargin{9};
  
  [alphaP, rhoP] =...
  wfnParamEstAlphaRho(FLaplace,TimeInput,NLag,alphamin,alphamax,rhomax,SearchSwitch,tolalpha,tolrho);
  
 otherwise
  error('wfnWeeksCoreAlphaRho','Incorrect number of inputs, 5 or 9 are required.');
end

%Laguerre polynomial coefficients calculation 
LaguerreCoef = wfncpuFFTLagCoefAlphaRho(FLaplace,NLag,alphaP,rhoP);

%Clenshaw algorithm to sum the polynomials
Invertf = wfnClenshawAlphaRho(NLag,TimeInput,alphaP,rhoP,LaguerreCoef);
 
%Estimate the absolute and relative error of the inverted function 
%Note: Reversing order of rho,alpha for this function is correct.
[AbsTotalError, AbsTruncateError, AbsRoundoffError] = wfncpuErrorEstAlphaRho(rhoP,alphaP,FLaplace,TimeInput,NLag,0); 

%Convert from the log base 10 value
AbsTotalError = 10^(AbsTotalError);
AbsTruncateError = 10^(AbsTruncateError);
AbsRoundoffError = 10^(AbsRoundoffError);

RelTotalError = abs(AbsTotalError)/abs(Invertf);

%This plots the trasformations of the contour (optional)
%fnSequenceTransforms(alphaP,rhoP,101)

end %function definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
