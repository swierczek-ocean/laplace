%myTestAdaptiveWeeks An adaptive integration based Weeks method
%   A Weeks method implementation that utilizes the internal MATLAB
%   adaptive integration code instead of a fast Fourier transform.
%   Note that the way points are known because, by construction,
%   the contour is the unit circle.
%   Patrick Kano, Moysey Brio
%   November 1, 2016
clear all;
clc;

format compact;
format long g;

TestCase = 2; %1 or 2

if(TestCase==1)
 %Example: F(s) = 1./s --> f(t) = 1.0
 Fofs = '(1./s)'
 Timevec = [1.0];
 alphaP = 0.000001;
 NLag = 32;
 ftruth = 1.0;
else   
 %Example: F(s) = '1./sqrt(s.^2+1)' --> f(t) = 'besselj(0,t)'
 Fofs = '1./sqrt(s.^2+1)'
 Timevec = 1.0;
 alphaP = 0.000001;
 NLag = 16;
 ftruth = besselj(0,1.0); %= 0.765197686557967  
end

%  Compare with
%   WeeksMethod('1./s',TimeInput,0.001,2)
%   WeeksMethod('1./sqrt(s.^2+1)',TimeInput,0.001,2)
%----------------------------%
%    sigma = (rho/2) - alpha;
%    b     = rho/2
%    alpha = b - sigma; 
%    rho   = 2b;

rhoP = 1.0; %Fixed!!!! Do not change
bP = rhoP/2.0;
sigmaP = bP - alphaP;

%Provided alpha (Just a a guess by me!)
tic;
InvertfA = wfnWeeksCoreAdaptiveIntegrate(@(s)eval(Fofs),Timevec,NLag,alphaP)
toc;

disp(['Invert A : ',num2str(InvertfA)]);

%Just test the search algorithm
RelErrorTol = 0.1
SearchSwitch = 5; 

%Optimized alpha internally using fminsearch
tic;
[InvertfB,alphaEst,rhoEst] = WeeksMethod(Fofs,Timevec,RelErrorTol,SearchSwitch)
toc;

%Diagnostics
disp(['Truth : Fixed Alpha : Searched Alpha ']);
[ftruth,InvertfA,InvertfB]

MeasAbsError = abs(InvertfB-ftruth)
MeasRelError = abs(InvertfB-ftruth)./abs(ftruth)