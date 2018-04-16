%Compare the Laguerre Expansion Coefficiets for the Sample Alpha,Rho
%Patrick Kano, Moysey Brio
%June 13, 2016
clear all;
close all;
clc;

%  A) Bessel function
  FLaplace = '1./sqrt(s.^2+1)'
%  f(t) = 'besselj(0,t)'
%
%  B) A-S_29.3.019
%  FLaplace = '1./(s.*(s.*s+9/25))'
%  f(t) = '(25/9)*(1-cos(3t/5))'
%
%  C) A-S_29.3.079
%  FLaplace = 'exp(1./s)./(s.^(1.5))'
%  f(t) = 'sinh(2*sqrt(t))/sqrt(pi)'
%
%  D) A-S_29.3.011 (challenging case for numerical inversion, not finite at t=0)
%  FLaplace = '(s+0.6).^(-9/11)'
%  f(t) = '(t.^(-2/11))./(gamma(9/11)*exp(0.6*t))'

%-------------------%
alphaP = 0.05;
rhoP = 1.0;
NLag = 25;


Fofs = @(s)eval(FLaplace)

tic;
AdaptLC =...
integral(@(w)wfnAdaptiveWeeksCore(w,Fofs,alphaP,rhoP,1,NLag),1,1,...
'ArrayValued',true,...
'Waypoints',[1i,-1,-1i],...
'AbsTol',0.001,...
'RelTol',0.001)
RunTimes.AdaptLC = toc;

tic;
ARfftLC = wfncpuFFTLagCoefAlphaRho(FLaplace, NLag, alphaP, rhoP)
RunTimes.ARfftLC = toc;

sigmaP = (rhoP/2.0) - alphaP;
bP = rhoP/2.0;
tic;
SPfftLC = wfncpuFFTLagCoefSigmab(FLaplace, NLag, sigmaP, bP)
RunTimes.SPfftLC = toc;

%not a vectorized
tic;
for nidx=1:NLag
 GaussKronrodLC(nidx) =....
 quadgk(@(w)wfnAdaptiveWeeksCore(w,Fofs,alphaP,rhoP,nidx,nidx),1,1,...
'Waypoints',[1i,-1,-1i],...
'AbsTol',0.001,...
'RelTol',0.001)
end %nidx
RunTimes.GKLC = toc;

tic;
TraciUJ=...
integral(@(theta)wfnTRACIcoreUJ(theta,Fofs,alphaP,rhoP,0,NLag-1),-pi,pi,...
'ArrayValued',true,...
'AbsTol',0.001,...
'RelTol',0.001);

TraciUH=...
integral(@(theta)wfnTRACIcoreUH(theta,Fofs,alphaP,rhoP,0,NLag-1),-pi,pi,...
'ArrayValued',true,...
'AbsTol',0.001,...
'RelTol',0.001);

TraciVJ=...
integral(@(theta)wfnTRACIcoreVJ(theta,Fofs,alphaP,rhoP,0,NLag-1),-pi,pi,...
'ArrayValued',true,...
'AbsTol',0.001,...
'RelTol',0.001);

TraciVH=...
integral(@(theta)wfnTRACIcoreVH(theta,Fofs,alphaP,rhoP,0,NLag-1),-pi,pi,...
'ArrayValued',true,...
'AbsTol',0.001,...
'RelTol',0.001);

 TraciLC = (TraciUJ + TraciUH) + 1i*(TraciVJ + TraciVH);
 
RunTimes.TraciLC = toc;

%-------------------------%
figure(101);
hold on;
plot(1:NLag,abs(AdaptLC),'-ko');
plot(1:NLag,abs(ARfftLC),'-rx');
plot(1:NLag,abs(SPfftLC),'-bd');
plot(1:NLag,abs(GaussKronrodLC),'-gh');
plot(1:NLag,abs(TraciLC),'-cp');
hold off;
legend('Adaptive','Alpha-Rho FFT','Sigma-b FFT','Gauss-Kronrod','TRACI');

figure(102);
subplot(1,2,1);
hold on;
plot(1:NLag,real(AdaptLC),'-ko');
plot(1:NLag,real(ARfftLC),'-rx');
plot(1:NLag,real(SPfftLC),'-bd');
plot(1:NLag,real(GaussKronrodLC),'-gh');
plot(1:NLag,real(TraciLC),'-cp');
hold off;

subplot(1,2,2);
hold on;
plot(1:NLag,imag(AdaptLC),'-ko');
plot(1:NLag,imag(ARfftLC),'-rx');
plot(1:NLag,imag(SPfftLC),'-bd');
plot(1:NLag,imag(GaussKronrodLC),'-gh');
plot(1:NLag,imag(TraciLC),'-cp');
hold off;
legend('Adaptive','Alpha-Rho FFT','Sigma-b FFT','Gauss-Kronrod','TRACI');

disp('RunTimes');
disp(RunTimes);

figure(103);
subplot(1,4,1);
plot(1:NLag,TraciUJ,'-ko');
legend('UJ');

subplot(1,4,2);
plot(1:NLag,TraciUH,'-ko');
legend('UH');

subplot(1,4,3);
plot(1:NLag,TraciVJ,'-ko');
legend('VJ');

subplot(1,4,4);
plot(1:NLag,TraciVH,'-ko');
legend('VH');
 