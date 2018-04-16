%Patrick Kano
%April 27, 2016
%Map the Bromwich contour back to the Riemann sphere
%    sigma = (rho/2) - alpha;
%    b     = rho/2
%    alpha = b - sigma; 
%    rho   = 2b;

%Correct but not used
%Rvec = sqrt([xvec.^2+yvec.^2+zvec.^2]);
%Thetavec = acos(zvec./Rvec);
%Phivec = atan(yvec./xvec);


clear all;
close all;
format compact;
format long g;


inFs = '(1./s).*(s-0.5 + 1)';

%bval = 10;
alphaval = 1;
%rho=b/2


%Ncs = 41;
%ColorMatrix = colormap(jet(Ncs));
%sigvec = linspace(-10,10,Ncs); 
Ncs = 1;
ColorMatrix = [0 0 1];
sigvec = 0.5;


figure(1);
hold on;
[xsph,ysph,zsph] = sphere(40);
Cmap=zeros(size(zsph));
surf(xsph,ysph,zsph,Cmap);
colormap('gray');
grid off;
shading interp;

figure(2);
hold on;
[xsph,ysph,zsph] = sphere(40);
Cmap=zeros(size(zsph));
surf(xsph,ysph,zsph,Cmap);
colormap('gray');
grid off;
shading interp;

tvec = linspace(-pi,pi,4096); %This is the 2D in W-plane theta

for sidx=1:length(sigvec)
sigma = sigvec(sidx)
bval = alphaval + sigma
alphaval=alphaval
%alphaval = bval-sigma %constant alpha test
rhoval = 2.0*bval

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

%Shift
Sone = Szero + alphaval;
svec = Sone;
denom = svec.*conj(svec)+1;
Xone = 2*real(svec)./denom;
Yone = 2*imag(svec)./denom;
Zone = (svec.*conj(svec)-1)./denom;
SphereThetaOne = acos(Zone); %r=1
SpherePhiOne = atan2(Yone,Xone);

%Invert
Stwo = 1./Sone;
svec = Stwo;
denom = svec.*conj(svec)+1;
Xtwo = 2*real(svec)./denom;
Ytwo = 2*imag(svec)./denom;
Ztwo = (svec.*conj(svec)-1)./denom;
SphereThetaTwo = acos(Ztwo); %r=1
SpherePhiTwo = atan2(Ytwo,Xtwo);

%Expand
Sthree = rhoval*Stwo;
svec = Sthree;
denom = svec.*conj(svec)+1;
Xthree = 2*real(svec)./denom;
Ythree = 2*imag(svec)./denom;
Zthree = (svec.*conj(svec)-1)./denom;
SphereThetaThree = acos(Zthree); %r=1
SpherePhiThree = atan2(Ythree,Xthree);

%Reflect
Sfour = -1*Sthree;
svec = Sfour;
denom = svec.*conj(svec)+1;
Xfour = 2*real(svec)./denom;
Yfour = 2*imag(svec)./denom;
Zfour = (svec.*conj(svec)-1)./denom;
SphereThetaFour = acos(Zfour); %r=1
SpherePhiFour = atan2(Yfour,Xfour);

%Shift up
Sfive = 1 + Sfour;
svec = Sfive;
denom = svec.*conj(svec)+1;
Xfive = 2*real(svec)./denom;
Yfive = 2*imag(svec)./denom;
Zfive = (svec.*conj(svec)-1)./denom;
SphereThetaFive = acos(Zfive); %r=1
SpherePhiFive = atan2(Yfive,Xfive);

svec = Wweeks;
denom = svec.*conj(svec)+1;
Xweeks = 2*real(svec)./denom;
Yweeks = 2*imag(svec)./denom;
Zweeks = (svec.*conj(svec)-1)./denom;
SphereThetaWeeks = acos(Zweeks); %r=1
SpherePhiWeeks = atan2(Yweeks,Xweeks);

%Sfinal should be the same as Weeks
disp('S--> W Transformations');
disp([Szero'; Sone'; Stwo'; Sthree'; Sfour'; Sfive'; Wweeks']);

%-------------------------------------------------------------------------%
figure(1);
plot3(Xzero,Yzero,Zzero,'-o','Color',ColorMatrix(sidx,1:3),'linewidth',1);

axis('square');
xlabel('x');
ylabel('y');
zlabel('z');

%-------------------%
figure(2);
hold on;
plot3(Xzero,Yzero,Zzero,'-ko','linewidth',2);
plot3(Xone,Yone,Zone,'-bo','linewidth',2);
plot3(Xtwo,Ytwo,Ztwo,'-go','linewidth',2);
plot3(Xthree,Ythree,Zthree,'-yo','linewidth',2);
plot3(Xfour,Yfour,Zfour,'-co','linewidth',2);
plot3(Xfive,Yfive,Zfive,'-mo','linewidth',2);
plot3(Xweeks,Yweeks,Zweeks,'-r+','linewidth',1);
hold off;
axis('square');
xlabel('x');
ylabel('y');
zlabel('z');
legend('Riemann Sphere',...
'S-Contour',...
'1: Translate by \alpha',...
'2: Inversion',...
'3: Dilation by \rho',...
'4: Rotation by \theta',...
'5: Translation by \beta',...
'W-Contour');

%-------------------------------%
figure(3);
subplot(1,6,1);
plot(real(Szero),imag(Szero),'-ko','linewidth',2);
axis square;
xlabel('S-Contour');

subplot(1,6,2);
plot(real(Sone),imag(Sone),'-bo','linewidth',2);
axis square;
xlabel('1: Translate by \alpha');

subplot(1,6,3);
plot(real(Stwo),imag(Stwo),'-go','linewidth',2);
axis square;
xlabel('2: Inversion');

subplot(1,6,4);
plot(real(Sthree),imag(Sthree),'-yo','linewidth',2);
axis square;
xlabel('3: Dilation by \rho');

subplot(1,6,5);
plot(real(Sfour),imag(Sfour),'-co','linewidth',2);
axis square;
xlabel('4: Rotation by \theta');

subplot(1,6,6);
hold on;
plot(real(Sfive),imag(Sfive),'-mo','linewidth',2);
plot(real(Wweeks),imag(Wweeks),'-r+','linewidth',1);
hold off;
axis square;
xlabel('5: Translation by \beta');
%'W-Contour');
drawnow;

%-----------------%
figure(10);
subplot(2,1,1);
hold on;
plot(real(Szero),imag(Szero),'-+','Color',ColorMatrix(sidx,1:3),'linewidth',0.1);
axis tight;

subplot(2,1,2);
hold on;
plot(real(Wweeks),imag(Wweeks),'-o','Color',ColorMatrix(sidx,1:3),'linewidth',1);
axis square;
drawnow;
axis tight;

figure(1);
title(['\alpha = ',num2str(alphaval),' \rho = ',num2str(rhoval), ' \sigma = ',num2str(sigma),' b = ',num2str(bval)]);

figure(10);
title(['\alpha = ',num2str(alphaval),' \rho = ',num2str(rhoval), ' \sigma = ',num2str(sigma),' b = ',num2str(bval)]);
pause(0.5);


%-----------------%
figure(11);
%Szero
%Xzero, Yzero, Zzero
%Fsp
subplot(1,2,1);
plot(tvec,abs(Fsp),'-ko');
ylabel('|F(s)|');

subplot(1,2,2);
plot(tvec,angle(Fsp),'-ko');
ylabel('phase(F(s))');

%-----------------%
figure(12);
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
  'Marker','none', ...
  'LineWidth',2);
view(3);
xlabel('real(W)');
ylabel('imag(W)');
zlabel('ln|F(s)|');

%-----------------%
figure(13);
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
  'Marker','none', ...
  'LineWidth',2);
view(3);
xlabel('Riemann Sphere X');
ylabel('Riemann Sphere Y');
ylabel('Riemann Sphere Z');
title('ln|F(s)|');

%-----------------%
figure(14);
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
xlabel('Riemann Sphere X');
ylabel('Riemann Sphere Y');
ylabel('Riemann Sphere Z');
title('angle|F(s)|');


%-----------------%
figure(15);
%Szero
%Xzero, Yzero, Zzero
%Fsp
Wz = log(abs(Fsp));
Wc = angle(Fsp);
h = surface(...
  'XData',[SphereThetaZero(:) SphereThetaZero(:)],...
  'YData',[SpherePhiZero(:) SpherePhiZero(:)],...
  'ZData',[Wz(:) Wz(:)],...
  'CData',[c(:) c(:)],...
  'FaceColor','none',...
  'EdgeColor','flat',...
  'Marker','none', ...
  'LineWidth',2);
view(3);
xlabel('Riemann Sphere \theta');
ylabel('Riemann Sphere \phi');
zlabel('ln|F(s)|');


end

disp('Complete');