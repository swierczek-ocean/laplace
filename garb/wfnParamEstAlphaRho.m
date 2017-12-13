function [alphaP, rhoP] = wfnParamEstAlphaRho(FLaplace,TimeInput,NLag,alphamin,alphamax,rhomax,SearchSwitch,tolalpha,tolrho)
%WFNPARAMESTALPHARHO Local and global searchs for optimal Weeks parameters sigma and b  
%  The function implements a local search based on the nested 1D fminbnd approach from J.A.C. Weideman
%  and a global search based on a uniform grid sampling in the (alpha,rho) plane.
% 
%  Use:
%  [alphaP, rhoP] = wfnParamEstAlphaRho(FLaplace,TimeInput,NLag,alphamin,alphamax,rhomax,SearchSwitch,tolalpha,tolrho)
%
%  Input:
%  FLaplace = a symbolic expression for the Laplace transform space function F(s)
%  TimeInput = time where f(t) is to be computed
%  NLag = number of terms in the Weeks expansion.
%  alphamin = minimum allowed value of alpha
%  alphamax = maximum allowed value of alpha
%  rhomax = maximum allowed value of rho
%  SearchSwitch = 0-CPU local fminbnd, 1-CPU global, 2-GPU/JACKET global 
%  tolalpha = tolerance on the alpha parameter
%  tolrho = tolerance on the rho parameter
%
%  Output:
%  alphaP = the optimized alpha parameter
%  rhoP = the optimized rho parameter
%
%  Author: 
%  Patrick Kano, Moysey Brio - 2016
%
%  Modification Date [M/D/Y]:
%  03/04/2016 - Version 1.0 based on previous sigma,b 
%  04/05/2016 - Updated with a fminsearch 2D
%  06/10/2016 - Removed option 5 and made it an adaptive integration

%   Options for the search switch are:
%   2 - CPU local fminsearch over alpha and rho
%   3 - CPU local fminbnd search over alpha and rho
%   4 - CPU global search over alpha and rho

DebugSwitch=1;
switch(SearchSwitch)
 case 2 %CPU local fminsearch
  ErrorFlag=0; %0-total error,1-truncation error only
  %ARproblem.object = '@(y)wfncpuErrorEstVecAR(y,FLaplace,TimeInput,NLag,ErrorFlag)';
  tmpoptions = optimset('TolX',0.0001*norm([tolrho,tolalpha]),'Display','iter');
 
  %sigma=0 --> alpha = rho/2
  [tmpAR,tmpErrorEst,tmpExitFlag,tmpOutput] =...
  fminsearch(@wfncpuErrorEstVecAR,[1,0],tmpoptions,FLaplace,TimeInput,NLag,ErrorFlag);
 
  rhoP = tmpAR(1)
  alphaP = tmpAR(2)
 
 case 3 %CPU local fminbnd
 
 alphaP = fminbnd('wfnNestedErrorAlphaRho',alphamin,alphamax,[0 tolalpha],FLaplace,TimeInput,NLag,rhomax,tolrho);
 rhoP = fminbnd('wfncpuErrorEstAlphaRho',0,rhomax,[0 tolrho],alphaP,FLaplace,TimeInput,NLag,0); %ErrorFlag=0->total error

 case 4 %global grid of alpha,rho
 
 Nalpha = floor((alphamax-alphamin)/tolalpha);
 Nrho = floor(rhomax/tolrho);
 
 Nprod = Nalpha*Nrho; 

 alphatemp = linspace(alphamin,alphamax,Nalpha+1); alphatemp = alphatemp(2:Nalpha+1); 
 rhotemp = linspace(0,rhomax,Nrho+1); rhotemp = rhotemp(2:Nrho+1);
  
 [alphamesh,rhomesh] = meshgrid(alphatemp,rhotemp);
 
 alphavec = alphamesh(:);
 rhovec = rhomesh(:);
 
 Errorvec = zeros(1,Nprod);
 
 for nidx=1:Nprod
  Errorvec(nidx) = wfncpuErrorEstAlphaRho(rhovec(nidx),alphavec(nidx),FLaplace,TimeInput,NLag,1);
 end

 [~, minidx] = min(Errorvec);
 alphaP = alphavec(minidx)
 rhoP = rhovec(minidx)
  
 %keyboard;
 if(DebugSwitch==1)
  figure(2001); 
  pcolor(alphamesh,rhomesh,reshape(Errorvec,Nrho,Nalpha));
  shading interp;
  xlabel('\alpha');
  ylabel('\rho'); 
  
  %Create a corresponding sigma,rho mesh
  figure(2002);
  subplot(1,2,1);
  pcolor(alphamesh,rhomesh,reshape(Errorvec,Nrho,Nalpha));
  shading flat;
  xlabel('\alpha');
  ylabel('\rho'); 
 
  sigmamesh = rhomesh/2.0-alphamesh;
  bmesh = rhomesh/2.0;
  subplot(1,2,2);
  pcolor(sigmamesh,bmesh,reshape(Errorvec,Nrho,Nalpha));
  shading flat;
  xlabel('\sigma');
  ylabel('b'); 
  drawnow;
 end
 
 otherwise 
  error('wfnParamEstAlphaRho:SearchSwitch','An incorrect choice for the search method has been made.');
end %switch

end %function definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
