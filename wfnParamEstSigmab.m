function [sigmaP, bP] = wfnParamEstSigmab(FLaplace,TimeInput,NLag,sigmin,sigmax,bmax,SearchSwitch,tols,tolb)
%WFNPARAMESTSIGMAB Local and global searchs for optimal Weeks parameters sigma and b  
%  The function implements a local search based on the nested 1D fminbnd approach from J.A.C. Weideman (method 2)
%  and a global search based on a uniform grid sampling in the (sigma,b) plane.
% 
%  Use:
%  [sigmaP, bP] = wfnParamEstSigmab(FLaplace,TimeInput,NLag,sigmin,sigmax,bmax,SearchSwitch,tols,tolb)
%
%  Input:
%  FLaplace = a symbolic expression for the Laplace transform space function F(s)
%  TimeInput = time where f(t) is to be computed
%  NLag = number of terms in the Weeks expansion.
%  sigmin = minimum allowed value of sigma
%  sigmax = maximum allowed value of sigma
%  bmax = maximum allowed value of b
%  SearchSwitch = 0-CPU local fminbnd, 1-CPU global
%  tols = tolerance on the s parameter
%  tolb = tolerance on the b parameter
%
%  Output:
%  sigmaP = the optimized sigma parameter
%  bP = the optimized b parameter
%
%  Author: 
%  Patrick Kano, Moysey Brio - 2011
%
%  Modification Date [M/D/Y]:
%  03/30/2011 - Version 1.0
%  03/04/2016 - Name change and Jacket GPU search removed

DebugSwitch=1;
switch SearchSwitch
 case 0 %CPU local fminbnd
 
  options = optimset('TolX',tols);
  sigmaP = fminbnd('wfnNestedErrorSigmab',sigmin,sigmax,options,FLaplace,TimeInput,NLag,bmax,tolb);
  options = optimset('TolX',tolb);
  bP = fminbnd('wfncpuErrorEstSigmab',0,bmax,options,sigmaP,FLaplace,TimeInput,NLag,0); %ErrorFlag=0->total error
  
     
%   sigmaP = fminbnd('wfnNestedErrorSigmab',sigmin,sigmax,[0 tols],FLaplace,TimeInput,NLag,bmax,tolb);
%   bP = fminbnd('wfncpuErrorEstSigmab',0,bmax,[0 tolb],sigmaP,FLaplace,TimeInput,NLag,0); %ErrorFlag=0->total error

 case 1 %CPU global

  Nsigma = floor((sigmax-sigmin)/tols);
  Nb = floor(bmax/tolb);
 
  Nprod = Nsigma*Nb; 
 
  sigmatemp = linspace(sigmin,sigmax,Nsigma+1); sigmatemp = sigmatemp(2:Nsigma+1); 
  btemp = linspace(0,bmax,Nb+1); btemp = btemp(2:Nb+1);
  
  [smesh,bmesh] = meshgrid(sigmatemp,btemp);
 
  svec = smesh(:);
  bvec = bmesh(:);
    
  Errorvec = zeros(1,Nprod);

  for nidx=1:Nprod
   Errorvec(nidx) = wfncpuErrorEstSigmab(bvec(nidx),svec(nidx),FLaplace,TimeInput,NLag,0);
  end

  [~, minidx] = min(Errorvec);
  sigmaP = svec(minidx);
  bP = bvec(minidx);

  %keyboard;
  if(DebugSwitch==1)
   figure(3001); 
   pcolor(smesh,bmesh,reshape(Errorvec,Nb,Nsigma))
   xlabel('\sigma');
   ylabel('b');
   
  %Create a corresponding sigma,rho mesh
  figure(3002);
  
  alphamesh = bmesh-smesh;
  rhomesh = 2.0*bmesh;
  subplot(1,2,1);
  pcolor(alphamesh,rhomesh,reshape(Errorvec,Nb,Nsigma));
  shading flat;
  xlabel('\alpha');
  ylabel('\rho'); 
 
  subplot(1,2,2);
  pcolor(smesh,bmesh,reshape(Errorvec,Nb,Nsigma));
  shading flat;
  xlabel('\sigma');
  ylabel('b'); 
  drawnow;
   
  end
  
  otherwise
   error('wfnParamEstSigmab:SearchSwitch','An incorrect choice for the search method has been made.');
 end %switch
end %function definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
