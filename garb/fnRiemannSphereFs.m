function fnRiemannSphereFs(inFs)
%FNRIEMANNSPHEREFS Create a color of F(s) on the Riemann Unit Sphere
%   Given Laplace space function F(s) as a string fnRiemannSphereFs(inFs),
%   this function creates stereographic
%   projections of the Riemann sphere on the complex plane for F(s)
%   In particular, it creates
%   a plot of log|F(s)| for s_{+}=(x+iy)/(1-z) in the usual 2D plane
%   and log|F(s)| for s_{-}=(x-iy)/(1+z) in the usual 2D plane
%   It also creates a 3D plot of |F(s)| on the unit sphere where
%   for the Northern hemisphere (x,y,z), z>0
%   s_{+} = s_{+}=(x+iy)/(1-z) = cot(\phi/2)e^{i\theta}
%   (x,y,z) = (2*Re(s), 2*Im(s), -1+conj(s)*s)/(1+s*conj(s))
%   and for the Southern hemisphere (x,y,z), z<0
%   s_{-}=(x-iy)/(1+z) = tan(phi/2)*exp^{-i\theta)
%   (x,y,z) = (2*Re(s), -2*Im(s), 1-conj(s)*s)/(1+s*conj(s))
%   As usual, elevation 0<=phi<=pi, azimuth 0<=\theta<2*pi
%
%   Patrick Kano, Moysey Brio
%   April 11, 2016
%   Example: fnRiemannSphereFs('1./s');
narginchk(1,1);

phivec = linspace(0,pi,1025); 
phivec = phivec(1:end-1);
thetavec = linspace(0,2*pi,1025);
thetavec = thetavec(1:end-1);
[phimesh,thetamesh] = meshgrid(phivec,thetavec);

sp = cot(0.5*phimesh).*exp(1i*thetamesh);
denomp = real(1+sp*conj(sp));
xp = (2*real(sp))./denomp;
yp = (-2*imag(sp))./denomp;
zp = real(sp.*conj(sp) - ones(size(sp)))./denomp;

s=sp;
Fsp = eval(inFs);

sn = tan(0.5*phimesh).*exp(-1i*thetamesh);
denomn = real(1+sn*conj(sn));
xn = (2*real(sn))./denomn;
yn = (2*imag(sn))./denomn;
zn = real(ones(size(sp))-sn.*conj(sn))./denomn;

s=sn;
Fsn = eval(inFs);

figure(1);
subplot(1,2,1);
surf(real(sp),imag(sp),log(abs(Fsp)));
axis square;
title('|F(s_{+})|');
shading interp;

subplot(1,2,2);
surf(real(sn),imag(sn),log(abs(Fsn)));
axis square;
title('|F(s_{-})|');
shading interp;

figure(2);
scatter3(xp(:),yp(:),zp(:),1+abs(Fsp(:)),angle(Fsp(:)));
title('F(s_{+})');
xlabel('x');
ylabel('y');
zlabel('z');
shading interp;
colorbar;
daspect([1 1 1]); 
axis equal;
drawnow;

figure(3);
%surf(xn,yn,zn,log(abs(Fsn)));
scatter3(xn(:),yn(:),zn(:),1+abs(Fsn(:)),angle(Fsn(:)));
title('F(s_{-})');
xlabel('x');
ylabel('y');
zlabel('z');
shading interp;
colorbar;
daspect([1 1 1]); 
axis equal; 


end %function definition