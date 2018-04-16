function fnSequenceTransforms(alpha,rho,Npoints)
%FNSEQUENCETRANSFORMS Create Image of the Sequence of Transforms to S
%   Create a sequece of images from S to W given by the Weeks method 
%   selection for the Moebius transformation
%   Patrick Kano
%   March 20, 2016

%% CONTOURS
%Sampling is uniform in W, but plot forward
Wvec = exp(1i*linspace(-pi,pi,Npoints));

% Bromwhich Contour
%Svec = (0.5*rho - alpha) + iuvec/2;
%uvec = 1i*rho*(Wvec+1)./(Wvec-1);

Svec = rho./(1-Wvec)-alpha;

%Shift
Wone = Svec + alpha;

%Invert
Wtwo = 1./Wone;

%Expand
Wthree = rho*Wtwo;

%Reflect
Wfour = -1*Wthree;

%Shift up
Wfinal = 1 + Wfour;

%% PLOTS
figure(1776);
clf;

subplot(1,6,1);
plot(real(Svec),imag(Svec),'-bx');
axis square;
xlabel('S-plane');

subplot(1,6,2);
plot(real(Wone),imag(Wone),'-rx');
axis square;
xlabel('Shift left by \alpha');

subplot(1,6,3);
plot(real(Wtwo),imag(Wtwo),'-gx');
axis square;
xlabel('Invert');

subplot(1,6,4);
plot(real(Wthree),imag(Wthree),'-mx');
axis square;
xlabel('Expand by \rho');

subplot(1,6,5);
plot(real(Wfour),imag(Wfour),'-cx');
axis square;
xlabel('Reflect across the origin');

subplot(1,6,6);
plot(real(Wfinal),imag(Wfinal),'-kx');
axis square;
xlabel('Shift right by 1');
drawnow;


end %function definition