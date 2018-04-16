%Create a full Riemann surface plot of F(s)
%Patrick O. Kano
%May 22, 2016

clear all;
close all;
clc;
inFs = '1./sqrt((s-1).*(s-2))'

% Selected values to try
SampleMethod = 0; %0-Sigma,b 1-Alpha,Rho

sigvec = linspace(0,5,6);

bvec = linspace(0,10,6);
bvec = bvec(2:end);

SampleMethod = 0; %0-Sigma,b 1-Alpha,Rho

alphavec = linspace(0,5,6);

rhovec = linspace(0,5,6);
rhovec = rhovec(2:end);

%% SETUP FIGURES
monitorpos = get(0, 'MonitorPositions');

figure(100);
set(gcf, 'Position',[0.7*monitorpos(3) 0.1*monitorpos(4) 0.3*monitorpos(3:4)])

figure(200);
set(gcf, 'Position',[0.4*monitorpos(3) 0.3*monitorpos(4) 0.3*monitorpos(3:4)])
subplot(1,2,1);
hold on;
subplot(1,2,2);
hold on;

figure(300);
set(gcf, 'Position',[0.1*monitorpos(3) 0.6*monitorpos(4) 0.3*monitorpos(3:4)])
hold on;

tvec = linspace(-pi,pi,1024); %This is the 2D in W-plane theta

%sigma,b
%for sidx=1:length(sigvec)
% for bidx=1:length(bvec)
% sigma = sigvec(sidx);
% bval = bvec(bidx);
% alphaval = bval-sigma; %constant alpha test
% rhoval = 2.0*bval;

%alpha,rho
for aidx=1:length(alphavec)
 for ridx=1:length(rhovec)
  alphaval = alphavec(aidx); 
  rhoval = rhovec(ridx);    
     
  sigma = 0.5*rhoval - alphaval;
  bval = 0.5*rhoval


if(rhoval==0),
 error('rho can not be zero');
end

%u from solving backwards from uniform w=e^{itheta} sampling to y
uvec = bval*sin(tvec)./(1-cos(tvec)); %finite sampling here yields the approximation
Szero = sigma + 1i*uvec;

s=Szero;
Fsp = eval(inFs);
clear s;

%Weeks figured out the exact Moebius transformation that will manipulate the
%Riemann sphere so that no matter which contour you choose, the sphereographic projected
%W-space contour from the manipulated sphere is ALWAYS the unit circle
%MY question: Is Weeks mapping unique?  Is it the only transformation of the contour
%that will always lead to W-projection on the unit circle?

Wweeks = (Szero - sigma - bval)./(Szero - sigma + bval);

svec = Szero;
denom = svec.*conj(svec)+1;
Xzero = 2*real(svec)./denom;
Yzero = 2*imag(svec)./denom;
Zzero = (svec.*conj(svec)-1)./denom;
SphereThetaZero = acos(Zzero); %r=1
SpherePhiZero = atan2(Yzero,Xzero);

%-----------------%
figure(100);
%Szero
%Xzero, Yzero, Zzero
%Fsp
x = real(Szero);
y = imag(Szero);
z = log(abs(Fsp));
c = angle(Fsp);
h = surface(...
  'XData',[x(:) x(:)],...
  'YData',[y(:) y(:)],...
  'ZData',[z(:) z(:)],...
  'CData',[c(:) c(:)],...
  'FaceColor','none',...
  'EdgeColor','flat',...
  'Marker','none', ...
  'LineWidth',2);
view(3);
xlabel('real(W)');
ylabel('imag(W)');
zlabel('ln|F(s)|');

%-----------------%
figure(200);
subplot(1,2,1);
hold on;
%Szero
%Xzero, Yzero, Zzero
%Fsp
Fspmag = log(abs(Fsp));
h = surface(...
  'XData',[Xzero(:) Xzero(:)],...
  'YData',[Yzero(:) Yzero(:)],...
  'ZData',[Zzero(:) Zzero(:)],...
  'CData',[Fspmag(:) Fspmag(:)],...
  'FaceColor','none',...
  'EdgeColor','flat',...
  'Marker','.', ...
  'LineWidth',1);
view(3);
xlabel('X');
ylabel('Y');
ylabel('Z');
title('ln|F(s)|');

%-----------------%
figure(200);
subplot(1,2,2);
hold on;
%Szero
%Xzero, Yzero, Zzero
%Fsp
Fspangle = angle(Fsp);
h = surface(...
  'XData',[Xzero(:) Xzero(:)],...
  'YData',[Yzero(:) Yzero(:)],...
  'ZData',[Zzero(:) Zzero(:)],...
  'CData',[Fspangle(:) Fspangle(:)],...
  'FaceColor','none',...
  'EdgeColor','flat',...
  'Marker','none', ...
  'LineWidth',2);
view(3);
xlabel('X');
ylabel('Y');
ylabel('Z');
title('angle|F(s)|');

%-----------------%
figure(300);
%Szero
%Xzero, Yzero, Zzero
%Fsp
x = real(Wweeks);
y = imag(Wweeks);
z = log(abs(Fsp));
c = angle(Fsp);
h = surface(...
  'XData',[x(:) x(:)],...
  'YData',[y(:) y(:)],...
  'ZData',[z(:) z(:)],...
  'CData',[c(:) c(:)],...
  'FaceColor','none',...
  'EdgeColor','flat',...
  'Marker','.', ...
  'LineWidth',1);
view(3);
xlabel('real(W)');
ylabel('imag(W)');
zlabel('ln|F(s)|');
end %bidx
end %sidx

disp('Complete');