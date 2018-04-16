function [alphaP, rhoP] = wfnParamEstAdaptiveAlpha(FLaplace,NLag,alphamin,alphamax,tolalpha)
%WFNPARAMESTALPHARHO Local and global searchs for optimal Weeks parameter alpha
%  The function implements a local search based on the nested 1D fminbnd approach 
%  from J.A.C. Weideman and a global search based on a uniform grid sampling in alpha.
% 
%  Use:
%  [alphaP, rhoP] = wfnParamEstAdaptiveAlpha(FLaplace,NLag,alphamin,alphamax,tolalpha)
%
%  Input:
%  FLaplace = a symbolic expression for the Laplace transform space function F(s)
%  NLag = number of terms in the Weeks expansion.
%  alphamin = minimum allowed value of alpha
%  alphamax = maximum allowed value of alpha
%  tolalpha = tolerance on the alpha parameter
%
%  Output:
%  alphaP = the optimized alpha parameter
%  rhoP = 1.0 (fixed)
%
%  Author: 
%  Patrick Kano, Moysey Brio - 2016
%
%  Modification Date [M/D/Y]:
%  06/10/2016 - Initial release

%  Note that the last argument to the estimation fuction is the 
%  ErrorFlag
%   0 = log base 10 of the total error
%   1 = log base 10 of the truncation error only
rhoP = 1.0;

 
%Brute force global search will work but is slow!

 Nalpha = 200; %floor((alphamax-alphamin)/tolalpha)
 
 alphavec = linspace(alphamin,alphamax,Nalpha+1);
 alphavec = alphavec(2:Nalpha+1); 
   
 Errorvec = zeros(1,Nalpha);

 for nidx=1:Nalpha, 
     nidx = nidx
  Errorvec(nidx) = wfnErrorEstAdaptiveAlpha(alphavec(nidx),FLaplace,NLag,0);
 end

 [~, minidx] = min(Errorvec);
 alphaP = alphavec(minidx);
 
 figure(2001);
 hold on;
 plot(alphavec,Errorvec,'-bx');
 plot(alphaP,Errorvec(minidx),'ko');
 hold off;
keyboard;

 

 ErrorFlag=0; %0-total error,1-truncation error onlyc
 tmpoptions = optimset('TolX',tolalpha,'Display','iter');
 %tmpoptions = optimset('TolFun',0.1,'Display','iter','MaxIter',25,'PlotFcns',@optimplotfval);
 
 %[alphaP,tmpErrorEst,tmpExitFlag,tmpOutput] =...
 %fminsearch(@wfnErrorEstAdaptiveAlpha,0.45,tmpoptions,FLaplace,NLag,ErrorFlag);
 
 %works too
 [alphaP,tmpErrorEst,tmpExitFlag,tmpOutput] =...
 fminbnd(@wfnErrorEstAdaptiveAlpha,alphamin,alphamax,tmpoptions,FLaplace,NLag,ErrorFlag)

 %This plots the trasformations of the contour (optional)
 %fnSequenceTransforms(alphaP,rhoP,101)
 
end %function definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
