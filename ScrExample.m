%SCREXAMPLE 
%  A script to test the Weeks method on F(s) expressions with known f(t):
%
%  A) Bessel function
%  F(s) = '1./sqrt(s.^2+1)'
%  f(t) = 'besselj(0,t)'
%  Timevec = 10:10:60;
%
%  B) A-S_29.3.019
%  F(s) = '1./(s.*(s.*s+9/25))'
%  f(t) = '(25/9)*(1-cos(3t/5))'
%  Timevec = 1:2:11;
%
%  C) A-S_29.3.079
%  F(s) = 'exp(1./s)./(s.^(1.5))'
%  f(t) = 'sinh(2*sqrt(t))/sqrt(pi)'
%  Timevec = 1:2:11;
%
%  D) A-S_29.3.011 (challenging case for numerical inversion, not finite at t=0)
%  F(s) = '(s+0.6).^(-9/11)'
%  f(t) = '(t.^(-2/11))./(gamma(9/11)*exp(0.6*t))'
%  Timevec = 0.5:0.5:3;
%
%  E) Dean Duffy Numerial Methods (equation 7.2.20)
%  F(s) = @(s)fnDuffyExample(s)
%  f(t) = 1 (Need this)
%  Timevec = 0.5:0.5:2.5;
%
%  The script implements a few different tests.
%  The first tests demonstrate how to call the core Weeks method function without the wrapper.
%  The additional tests are of the code with the wrapper.
%
%  A-1) Core Weeks Method: without a search (sigma=1,b=1) 
%  A-2) Core Weeks Method: sigma-b CPU Local Search  
%  A-3) Core Weeks Method: sigma-b CPU Global Search 
%  A-4) Core Weeks Method: alpha-rho Search
%  A-5) Wrapped Weeks Method: sigma-b CPU Local Search
%  A-6) Wrapped Weeks Method: alpha-rho Global Search
%  A-7) Clenshaw Algorithm/Wrapped Weeks Method: sigma-b CPU Local Search
%  A-8) Wrapped Weeks Method: Adaptive Integrative alpha Fminsearch
%
%  B-1) Wrapped Weeks Method: sigma-b CPU Local Search
%  B-2) Wrapped Weeks Method: alpha-rho Search
%  B-3) Wrapped Weeks Method: Adaptive Integrative alpha Fminsearch
%
%  C-1) Wrapped Weeks Method: sigma-b CPU Local Search
%  C-2) Wrapped Weeks Method: alpha-rho Search
%  C-3) Wrapped Weeks Method: Adaptive Integrative alpha Fminsearch
%
%  D-1) Wrapped Weeks Method: sigma-b CPU Local Search
%  D-2) Wrapped Weeks Method: alpha-rho Search
%  D-3) Wrapped Weeks Method: Adaptive Integrative alpha Fminsearch
%
%  E-1) Wrapped Weeks Method: alpha-rho fminsearch Search
%  E-2) Wrapped Weeks Method: alpha-rho local fminbnd
%  E-3) Wrapped Weeks Method: global alpha-rho search 
%  E-4) Wrapped Weeks Method: Adaptive Integrative alpha Fminsearch
%
%  If it returns a fail value, this script will only compute the CPU examples.
%
%  References:
%  The numbers of the examples are those found in:
%  1) "Inversion of noise-free Laplace transforms: 
%  Towards a standardized set of test problems", P. Valko and S. Vadja
%  Inverse Problems in Engineering, volume 10, no. 5, pp. 467-483 (2002).  
%  www.pe.tamu.edu/valko/public%5Fhtml/Nil/Tests/index.htm
%  2) "Algorithms for Parameter Selection in the Weeks Method for Inverting 
%  the Laplace Transform", A. Weideman, SIAM J. on Scientific Computing, 
%  volume 21, pp. 111-128 (1999). 
%  http://dip.sun.ac.za/~weideman/research/weeks.html
%
%  Author: 
%  Patrick Kano, Moysey Brio - 2011, 2016
%
%  Modification Date [M/D/Y]:
%  03/30/2011 - Version 1.0
%  04/08/2016 - Version 2.0 - Removed Jacket tests
%  06/12/2016 - Version 2.1 - Added adapative Weeks method tests

clear all;
close all;
format compact;
format long g;

%1=0n,0=Off
SetTest.A1 = 0; %  A-1) Core Weeks Method: without a search (sigma=1,b=1) 
SetTest.A2 = 0; %  A-2) Core Weeks Method: sigma-b CPU Local Search  
SetTest.A3 = 0; %  A-3) Core Weeks Method: sigma-b CPU Global Search 
SetTest.A4 = 0; %  A-4) Core Weeks Method: alpha-rho Search
SetTest.A5 = 0; %  A-5) Wrapped Weeks Method: sigma-b CPU Local Search
SetTest.A6 = 0; %  A-6) Wrapped Weeks Method: alpha-rho Global Search
SetTest.A7 = 0; %  A-7) Clenshaw Algorithm/Wrapped Weeks Method: sigma-b CPU Local Search
SetTest.A8 = 0; %  A-8) Wrapped Weeks Method: Adaptive Integrative alpha Fminsearch
%
SetTest.B1 = 0; %  B-1) Wrapped Weeks Method: sigma-b CPU Local Search
SetTest.B2 = 0; %  B-2) Wrapped Weeks Method: alpha-rho Search
SetTest.B3 = 0; %  B-3) Wrapped Weeks Method: Adaptive Integrative alpha Fminsearch
%
SetTest.C1 = 0; %  C-1) Wrapped Weeks Method: sigma-b CPU Local Search
SetTest.C2 = 0; %  C-2) Wrapped Weeks Method: alpha-rho Search
SetTest.C3 = 1; %  C-3) Wrapped Weeks Method: Adaptive Integrative alpha Fminsearch
%
SetTest.D1 = 0; %  D-1) Wrapped Weeks Method: sigma-b CPU Local Search
SetTest.D2 = 0; %  D-2) Wrapped Weeks Method: alpha-rho Search
SetTest.D3 = 0; %  D-3) Wrapped Weeks Method: Adaptive Integrative alpha Fminsearch
%
SetTest.E1 = 0; %  E-1) Wrapped Weeks Method: alpha-rho fminsearch Search
SetTest.E2 = 0; %  E-2) Wrapped Weeks Method: alpha-rho fminbnd Search
SetTest.E3 = 0; %  E-3) Wrapped Weeks Method: alpha-rho global Search
SetTest.E4 = 0; %  E-4) Wrapped Weeks Method: Adaptive Integrative alpha Fminsearch
%%%%%%%%%%%%%%%
if(SetTest.A1==1)
disp('------------------------------------------------------------');
disp('Test A-1');
disp('Core Weeks Method: without a search (sigma=1,b=1)');

Fs = '1./sqrt(s.^2+1)'
ft = 'besselj(0,t)'
Timevec = 10:10:60;
sigmaPinput = ones(1,6);
bPinput = ones(1,6);

fexact = feval(inline(ft,'t'),Timevec);

Nvec = [8 16 32 64 128];

for nidx=1:length(Nvec)
 NLag = Nvec(nidx);
 
 for tidx=1:length(Timevec)
  tic;
  [fapprox(nidx,tidx),sigmaP(nidx,tidx),bP(nidx,tidx)] =...
  wfnWeeksCoreSigmab(Fs,Timevec(tidx),NLag,sigmaPinput(tidx),bPinput(tidx));
  Runtime(nidx,tidx) = toc;
 end %tidx
 
 MeaAbsError(nidx,:) = 100*abs(fexact-fapprox(nidx,:));    %percent
 MeaRelError(nidx,:) = MeaAbsError(nidx,:)./abs(fexact); %percent
end

disp('Time: ');  		disp(Timevec);
disp('Exact: '); 		disp(fexact); 
disp('Approximate: ');      	disp(real(fapprox)); 
disp('Absolute Error %: '); 	disp(MeaAbsError); 
disp('Relative Error %: '); 	disp(MeaRelError); 
disp('Run Time :');         	disp(Runtime);
disp('sigma: '); 		disp(sigmaP);
disp('b:');      		disp(bP);

figure(1);
title('A-1');
hold on;
plot(Timevec,fexact,'-ko','LineWidth',2);
plot(Timevec,real(fapprox(1,:)),'-rx');
plot(Timevec,real(fapprox(2,:)),'-gx');
plot(Timevec,real(fapprox(3,:)),'-bx');
plot(Timevec,real(fapprox(4,:)),'-cx');
plot(Timevec,real(fapprox(5,:)),'-mx');
hold off;
axis([0 60 -0.5 0.2]);
xlabel('Time [sec]');
legend('Exact','N=8','N=16','N=32','N=64','N=128');
drawnow;

saveas(gcf,[pwd,filesep,'EXAMPLE_FIGURES',filesep,'A_1.jpg']);
clearvars -except SetTest;
end 

%%%%%%%%%%%%%%%
if(SetTest.A2==1)
disp('------------------------------------------------------------');
disp('Test A-2');
disp('Core Weeks Method: CPU Local Search');

Fs = '1./sqrt(s.^2+1)'
ft = 'besselj(0,t)'
Timevec = 10:10:60;
sigmamin = 0;
sigmamax = 20;
bmax = 20;
SearchSwitch = 0; %CPU-local
tols = 0.2;
tolb = 0.2;

fexact = feval(inline(ft,'t'),Timevec);

Nvec = [8 16 32 64 128];

for nidx=1:length(Nvec)
 NLag = Nvec(nidx);
 
 for tidx=1:length(Timevec)
  tic;
  [fapprox(nidx,tidx),sigmaP(nidx,tidx),bP(nidx,tidx)] =...
  wfnWeeksCoreSigmab(Fs,Timevec(tidx),NLag,sigmamin,sigmamax,bmax,SearchSwitch,tols,tolb);
  Runtime(nidx,tidx) = toc;
 end %tidx

 MeaAbsError(nidx,:) = 100*abs(fexact-fapprox(nidx,:));    %percent
 MeaRelError(nidx,:) = MeaAbsError(nidx,:)./abs(fexact); %percent
end

disp('Time: ');  		disp(Timevec);
disp('Exact: '); 		disp(fexact); 
disp('Approximate: ');      	disp(real(fapprox)); 
disp('Absolute Error %: '); 	disp(MeaAbsError); 
disp('Relative Error %: '); 	disp(MeaRelError); 
disp('Run Time :');         	disp(Runtime);
disp('sigma: '); 		disp(sigmaP);
disp('b:');      		disp(bP);

figure(2);
title('A-2');
hold on;
plot(Timevec,fexact,'-ko','LineWidth',2);
plot(Timevec,real(fapprox(1,:)),'-rx');
plot(Timevec,real(fapprox(2,:)),'-gx');
plot(Timevec,real(fapprox(3,:)),'-bx');
plot(Timevec,real(fapprox(4,:)),'-cx');
plot(Timevec,real(fapprox(5,:)),'-mx');
hold off;
xlabel('Time [sec]');
legend('Exact','N=8','N=16','N=32','N=64','N=128','Location','southeast');
drawnow;

saveas(gcf,[pwd,filesep,'EXAMPLE_FIGURES',filesep,'A_2.jpg']);
clearvars -except SetTest;
end

%%%%%%%%%%%%%%%
if(SetTest.A3==1)
disp('------------------------------------------------------------');
disp('Test A-3');
disp('Core Weeks Method: CPU Global Local (slow)');

Fs = '1./sqrt(s.^2+1)'
ft = 'besselj(0,t)'
Timevec = 10:10:60;
sigmamin = 0;
sigmamax = 20;
bmax = 20;
SearchSwitch = 1; %CPU-global
tols = 0.2;
tolb = 0.2;

fexact = feval(inline(ft,'t'),Timevec);

Nvec = [8 16 32 64 128];

for nidx=1:length(Nvec)
 NLag = Nvec(nidx);
 
 for tidx=1:length(Timevec)
  tic;
  [fapprox(nidx,tidx),sigmaP(nidx,tidx),bP(nidx,tidx)] =...
  wfnWeeksCoreSigmab(Fs,Timevec(tidx),NLag,sigmamin,sigmamax,bmax,SearchSwitch,tols,tolb);
  Runtime(nidx,tidx) = toc;
 end %tidx

 MeaAbsError(nidx,:) = 100*abs(fexact-fapprox(nidx,:));    %percent
 MeaRelError(nidx,:) = MeaAbsError(nidx,:)./abs(fexact); %percent
end

disp('Time: ');  		disp(Timevec);
disp('Exact: '); 		disp(fexact); 
disp('Approximate: ');      	disp(real(fapprox)); 
disp('Absolute Error %: '); 	disp(MeaAbsError); 
disp('Relative Error %: '); 	disp(MeaRelError); 
disp('Run Time :');         	disp(Runtime);
disp('sigma: '); 		disp(sigmaP);
disp('b:');      		disp(bP);

figure(3);
title('A-3');
hold on;
plot(Timevec,fexact,'-ko','LineWidth',2);
plot(Timevec,real(fapprox(1,:)),'-rx');
plot(Timevec,real(fapprox(2,:)),'-gx');
plot(Timevec,real(fapprox(3,:)),'-bx');
plot(Timevec,real(fapprox(4,:)),'-cx');
plot(Timevec,real(fapprox(5,:)),'-mx');
hold off;
xlabel('Time [sec]');
legend('Exact','N=8','N=16','N=32','N=64','N=128');
drawnow;

saveas(gcf,[pwd,filesep,'EXAMPLE_FIGURES',filesep,'A_3.jpg']);
clearvars -except SetTest;
end

%%%%%%%%%%%%%%%
if(SetTest.A4==1)
disp('------------------------------------------------------------');
disp('Test A-4');
disp('Core Weeks Method: Alpha-Rho Search (fast)');

Fs = '1./sqrt(s.^2+1)'
ft = 'besselj(0,t)'
Timevec = 10:10:60;
alphamin = NaN;
alphamax = NaN;
rhomax = NaN;
SearchSwitch = 2; %fminsearch
tolalpha = 0.02;
tolrho = 0.02;

fexact = feval(inline(ft,'t'),Timevec);

Nvec = [8 16 32 64 128];

for nidx=1:length(Nvec)
 NLag = Nvec(nidx);
 
 for tidx=1:length(Timevec)
  tic;
  [fapprox(nidx,tidx),alphaP(nidx,tidx),rhoP(nidx,tidx)] =...
  wfnWeeksCoreAlphaRho(Fs,Timevec(tidx),NLag,alphamin,alphamax,rhomax,SearchSwitch,tolalpha,tolrho);
  Runtime(nidx,tidx) = toc;
 end %tidx

 MeaAbsError(nidx,:) = 100*abs(fexact-fapprox(nidx,:));    %percent
 MeaRelError(nidx,:) = MeaAbsError(nidx,:)./abs(fexact); %percent
end

disp('Time: ');  		disp(Timevec);
disp('Exact: '); 		disp(fexact); 
disp('Approximate: ');      	disp(real(fapprox)); 
disp('Absolute Error %: '); 	disp(MeaAbsError); 
disp('Relative Error %: '); 	disp(MeaRelError); 
disp('Run Time :');         	disp(Runtime);
disp('alpha: '); 		disp(alphaP);
disp('rho:');      		disp(rhoP);

figure(4);
title('A-4');
hold on;
plot(Timevec,fexact,'-ko','LineWidth',2);
plot(Timevec,real(fapprox(1,:)),'-rx');
plot(Timevec,real(fapprox(2,:)),'-gx');
plot(Timevec,real(fapprox(3,:)),'-bx');
plot(Timevec,real(fapprox(4,:)),'-cx');
plot(Timevec,real(fapprox(5,:)),'-mx');
hold off;
xlabel('Time [sec]');
legend('Exact','N=8','N=16','N=32','N=64','N=128');
drawnow;

saveas(gcf,[pwd,filesep,'EXAMPLE_FIGURES',filesep,'A_4.jpg']);
clearvars -except SetTest;
end

%%%%%%%%%%%%%%%
if(SetTest.A5==1)
disp('------------------------------------------------------------');
disp('Test A-5');
disp('Wrapped Weeks Method: CPU Local Search');

Fs = '1./sqrt(s.^2+1)'
ft = 'besselj(0,t)'
Timevec = 10:10:60;
RelErrorTol = [0.01 0.01 0.01 0.01 0.01 0.01]; %1 percent relative error

fexact = feval(inline(ft,'t'),Timevec);

tic;
[fapprox,sigmaP,bP,EstRelError,EstAbsError,EstTruncateError,EstRoundoffError,ToleranceFlag,LaguerreCoef] =...
WeeksMethod(Fs,Timevec,RelErrorTol);
Runtime = toc;

MeaAbsError = 100*abs(fexact-fapprox);    %percent
MeaRelError = MeaAbsError./abs(fexact); %percent

disp('Time: ');  			disp(Timevec);
disp('Exact: '); 			disp(fexact); 
disp('Approximate: ');      		disp(real(fapprox)); 
disp('Measured Absolute Error %: ');  	disp(MeaAbsError); 
disp('Estimated Absolute Error %: '); 	disp(EstAbsError*100);
disp('Measured Relative Error %: ');  	disp(MeaRelError); 
disp('Estimated Relative Error %: '); 	disp(EstRelError*100);
disp('Run Time :');        		disp(Runtime');
disp('sigma: '); 			disp(sigmaP);
disp('b:');      			disp(bP);
disp('Tolerance Flag: ');		disp(ToleranceFlag);

figure(5);
title('A-5');
hold on;
plot(Timevec,fexact,'-ko','LineWidth',2);
plot(Timevec,real(fapprox),'-rx');
hold off;
xlabel('Time [sec]');
legend('Exact','Approximate');
drawnow;

saveas(gcf,[pwd,filesep,'EXAMPLE_FIGURES',filesep,'A_5.jpg']);
clearvars -except SetTest;
end

%%%%%%%%%%%%%%%
if(SetTest.A6==1)
disp('------------------------------------------------------------');
disp('Test A-6');
disp('Wrapped Weeks Method: Alpha-Rho Search');

Fs = '1./sqrt(s.^2+1)'
ft = 'besselj(0,t)'
Timevec = 10:10:60;
RelErrorTol = [0.01 0.01 0.01 0.01 0.01 0.01]; %1 percent relative error

fexact = feval(inline(ft,'t'),Timevec);

tic;
[fapprox,alphaP,rhoP,EstRelError,EstAbsError,EstTruncateError,EstRoundoffError,ToleranceFlag,LaguerreCoef] =...
WeeksMethod(Fs,Timevec,RelErrorTol,4);
Runtime = toc;

MeaAbsError = 100*abs(fexact-fapprox);    %percent
MeaRelError = MeaAbsError./abs(fexact); %percent

disp('Time: ');  			disp(Timevec);
disp('Exact: '); 			disp(fexact); 
disp('Approximate: ');      		disp(real(fapprox)); 
disp('Measured Absolute Error %: ');  	disp(MeaAbsError); 
disp('Estimated Absolute Error %: '); 	disp(EstAbsError*100);
disp('Measured Relative Error %: ');  	disp(MeaRelError); 
disp('Estimated Relative Error %: '); 	disp(EstRelError*100);
disp('Run Time :');        		disp(Runtime');
disp('alpha: '); 			disp(alphaP);
disp('rho:');      			disp(rhoP);
disp('Tolerance Flag: ');		disp(ToleranceFlag);

figure(6);
title('A-6');
hold on;
plot(Timevec,fexact,'-ko','LineWidth',2);
plot(Timevec,real(fapprox),'-rx');
hold off;
xlabel('Time [sec]');
legend('Exact','Approximate');
drawnow;

saveas(gcf,[pwd,filesep,'EXAMPLE_FIGURES',filesep,'A_6.jpg']);
clearvars -except SetTest;
end

%%%%%%%%%%%%%%%
if(SetTest.A7==1)
disp('------------------------------------------------------------');
disp('Test A-7');
disp('Clenshaw Algorithm/Wrapped Weeks Method: CPU Local Search');

%This purpose of this example is to show that it is in principle possible to 
%use the Laguerre coefficients directly from the inversion at a single time
%to rapidly compute a coarser estimate of f(t) at multiple times.
%An alternative is to iterpolate (sigma,b) from a coarse time sampling to a finer time sampling
%and then evaluate the core Weeks Function at these values.
%The tool returns sufficient information to the user to investigate multiple options.

tic;
Fs = '1./sqrt(s.^2+1)'
ft = 'besselj(0,t)'
Timevec = 60;

[~,sigmaP,bP,~,~,~,~,~,LaguerreCoef] = WeeksMethod(Fs,Timevec);

Timeinterp = linspace(0,60,1000)'; %covert to a column vector

%Extract the vector of coefficients from the cell array
vecLaguerreCoef = LaguerreCoef{1}';
NLag = length(vecLaguerreCoef);

fapprox = wfnClenshawSigmab(NLag,Timeinterp,sigmaP,bP,vecLaguerreCoef);

%Compare with the analytic result
fexact = feval(inline(ft,'t'),Timeinterp);
MeaAbsError = 100*abs(fexact-fapprox);    %percent
MeaRelError = MeaAbsError./abs(fexact); %percent

Runtime = toc

figure(7);
title('A-7: Solution');
hold on;
plot(Timeinterp,fexact,'-ko','LineWidth',2);
plot(Timeinterp,real(fapprox),'-rx');
hold off;
xlabel('Time [sec]');
legend('Exact','Approximate');
drawnow;
saveas(gcf,[pwd,filesep,'EXAMPLE_FIGURES',filesep,'A_7_soln.jpg']);

figure(8);
plot(Timeinterp,log10(MeaRelError),'-ko');
title('A-7: Relative Error');
xlabel('Time [sec]');
ylabel('log10(Relative Error)');
drawnow;
saveas(gcf,[pwd,filesep,'EXAMPLE_FIGURES',filesep,'A_7_error.jpg']);

clearvars -except SetTest;
end

%%%%%%%%%%%%%%%
if(SetTest.A8==1)
disp('------------------------------------------------------------');
disp('Test A-8');
disp('Wrapped Weeks Method: Adaptive Integrative alpha Fminsearch');

Fs = '1./sqrt(s.^2+1)'
ft = 'besselj(0,t)'
Timevec = 10:10:60;
RelErrorTol = [0.01 0.01 0.01 0.01 0.01 0.01]; %1 percent relative error

fexact = feval(inline(ft,'t'),Timevec);

tic;
[fapprox,alphaP,rhoP,EstRelError,EstAbsError,EstTruncateError,EstRoundoffError,ToleranceFlag,LaguerreCoef] =...
WeeksMethod(Fs,Timevec,RelErrorTol,5);
Runtime = toc;

MeaAbsError = 100*abs(fexact-fapprox);    %percent
MeaRelError = MeaAbsError./abs(fexact); %percent

disp('Time: ');                         disp(Timevec);
disp('Exact: ');                        disp(fexact); 
disp('Approximate: ');                  disp(real(fapprox)); 
disp('Measured Absolute Error %: ');  	disp(MeaAbsError); 
disp('Estimated Absolute Error %: '); 	disp(EstAbsError*100);
disp('Measured Relative Error %: ');  	disp(MeaRelError); 
disp('Estimated Relative Error %: '); 	disp(EstRelError*100);
disp('Run Time :');                     disp(Runtime');
disp('alpha: ');                        disp(alphaP);
disp('rho:');                           disp(rhoP);
disp('Tolerance Flag: ');               disp(ToleranceFlag);

figure(9);
title('A-8');
hold on;
plot(Timevec,fexact,'-ko','LineWidth',2);
plot(Timevec,real(fapprox),'-rx');
hold off;
xlabel('Time [sec]');
legend('Exact','Approximate');
drawnow;

saveas(gcf,[pwd,filesep,'EXAMPLE_FIGURES',filesep,'A_8.jpg']);
clearvars -except SetTest;
end

%%%%%%%%%%%%%%%
if(SetTest.B1==1)
%B) A-S_29.3.019
disp('------------------------------------------------------------');
disp('Test B-1');
disp('Wrapped Weeks Method: CPU Local Search');

Fs = '1./(s.*(s.*s+9/25))'
ft = '(25/9)*(1-cos(3*t/5))'
Timevec = 1:2:11;
RelErrorTol = [0.01 0.01 0.01 0.01 0.01 0.01]; %1 percent relative error

fexact = feval(inline(ft,'t'),Timevec);

tic;
[fapprox,sigmaP,bP,EstRelError,EstAbsError,~,~,ToleranceFlag] = WeeksMethod(Fs,Timevec,RelErrorTol);
Runtime = toc;

MeaAbsError = 100*abs(fexact-fapprox);    %percent
MeaRelError = MeaAbsError./abs(fexact); %percent

disp('Time: ');  			disp(Timevec);
disp('Exact: '); 			disp(fexact); 
disp('Approximate: ');      		disp(real(fapprox)); 
disp('Measured Absolute Error %: ');  	disp(MeaAbsError); 
disp('Estimated Absolute Error %: '); 	disp(EstAbsError*100);
disp('Measured Relative Error %: ');  	disp(MeaRelError); 
disp('Estimated Relative Error %: '); 	disp(EstRelError*100);
disp('Run Time :');        		disp(Runtime');
disp('sigma: '); 			disp(sigmaP);
disp('b:');      			disp(bP);
disp('Tolerance Flag: ');		disp(ToleranceFlag);

figure(10);
title('B-1');
hold on;
plot(Timevec,fexact,'-ko','LineWidth',2);
plot(Timevec,real(fapprox),'-rx');
hold off;
xlabel('Time [sec]');
legend('Exact','Approximate');
drawnow;

saveas(gcf,[pwd,filesep,'EXAMPLE_FIGURES',filesep,'B_1.jpg']);

clearvars -except SetTest;
end

%%%%%%%%%%%%%%%
if(SetTest.B2==1)
%B) A-S_29.3.019
disp('------------------------------------------------------------');
disp('Test B-2');
disp('Wrapped Weeks Method: Alpha-Rho Global Search');

Fs = '1./(s.*(s.*s+9/25))'
ft = '(25/9)*(1-cos(3*t/5))'
Timevec = 1:2:11;
RelErrorTol = [0.01 0.01 0.01 0.01 0.01 0.01]; %1 percent relative error

fexact = feval(inline(ft,'t'),Timevec);

tic;
[fapprox,alphaP,rhoP,EstRelError,EstAbsError,~,~,ToleranceFlag] = WeeksMethod(Fs,Timevec,RelErrorTol,2);
Runtime = toc;

MeaAbsError = 100*abs(fexact-fapprox);    %percent
MeaRelError = MeaAbsError./abs(fexact); %percent

disp('Time: ');  			disp(Timevec);
disp('Exact: '); 			disp(fexact); 
disp('Approximate: ');      		disp(real(fapprox)); 
disp('Measured Absolute Error %: ');  	disp(MeaAbsError); 
disp('Estimated Absolute Error %: '); 	disp(EstAbsError*100);
disp('Measured Relative Error %: ');  	disp(MeaRelError); 
disp('Estimated Relative Error %: '); 	disp(EstRelError*100);
disp('Run Time :');        		disp(Runtime');
disp('alpha: '); 			disp(alphaP);
disp('rho:');      			disp(rhoP);
disp('Tolerance Flag: ');		disp(ToleranceFlag);

figure(11);
title('B-2');
hold on;
plot(Timevec,fexact,'-ko','LineWidth',2);
plot(Timevec,real(fapprox),'-rx');
hold off;
xlabel('Time [sec]');
legend('Exact','Approximate');
drawnow;

saveas(gcf,[pwd,filesep,'EXAMPLE_FIGURES',filesep,'B_2.jpg']);

clearvars -except SetTest;
end

if(SetTest.B3==1)
%B) A-S_29.3.019
disp('------------------------------------------------------------');
disp('Test B-3');
disp('Wrapped Weeks Method: Adaptive Integrative alpha Fminsearch');

Fs = '1./(s.*(s.*s+9/25))'
ft = '(25/9)*(1-cos(3*t/5))'
Timevec = 1:2:11;
RelErrorTol = [0.01 0.01 0.01 0.01 0.01 0.01]; %1 percent relative error

fexact = feval(inline(ft,'t'),Timevec);

tic;
[fapprox,alphaP,rhoP,EstRelError,EstAbsError,~,~,ToleranceFlag] = WeeksMethod(Fs,Timevec,RelErrorTol,5);
Runtime = toc;

MeaAbsError = 100*abs(fexact-fapprox);    %percent
MeaRelError = MeaAbsError./abs(fexact); %percent

disp('Time: ');  			disp(Timevec);
disp('Exact: '); 			disp(fexact); 
disp('Approximate: ');      		disp(real(fapprox)); 
disp('Measured Absolute Error %: ');  	disp(MeaAbsError); 
disp('Estimated Absolute Error %: '); 	disp(EstAbsError*100);
disp('Measured Relative Error %: ');  	disp(MeaRelError); 
disp('Estimated Relative Error %: '); 	disp(EstRelError*100);
disp('Run Time :');        		disp(Runtime');
disp('alpha: '); 			disp(alphaP);
disp('rho:');      			disp(rhoP);
disp('Tolerance Flag: ');		disp(ToleranceFlag);

figure(12);
title('B-3');
hold on;
plot(Timevec,fexact,'-ko','LineWidth',2);
plot(Timevec,real(fapprox),'-rx');
hold off;
xlabel('Time [sec]');
legend('Exact','Approximate');
drawnow;

saveas(gcf,[pwd,filesep,'EXAMPLE_FIGURES',filesep,'B_3.jpg']);

clearvars -except SetTest;
end

%%%%%%%%%%%%%%%
if(SetTest.C1==1)
% C) A-S_29.3.079
disp('------------------------------------------------------------');
disp('Test C-1');
disp('Wrapped Weeks Method: CPU Local Search');

Fs = 'exp(1./s)./(s.^(1.5))';
ft = 'sinh(2*sqrt(t))/sqrt(pi)'
Timevec = 1:2:11;
RelErrorTol = [0.01 0.01 0.01 0.01 0.01 0.01]; %1 percent relative error

fexact = feval(inline(ft,'t'),Timevec);

tic;
[fapprox,sigmaP,bP,EstRelError,EstAbsError,~,~,ToleranceFlag] = WeeksMethod(Fs,Timevec,RelErrorTol);
Runtime = toc;

MeaAbsError = 100*abs(fexact-fapprox);    %percent
MeaRelError = MeaAbsError./abs(fexact); %percent

disp('Time: ');  			disp(Timevec);
disp('Exact: '); 			disp(fexact); 
disp('Approximate: ');      		disp(real(fapprox)); 
disp('Measured Absolute Error %: ');  	disp(MeaAbsError); 
disp('Estimated Absolute Error %: '); 	disp(EstAbsError*100);
disp('Measured Relative Error %: ');  	disp(MeaRelError); 
disp('Estimated Relative Error %: '); 	disp(EstRelError*100);
disp('Run Time :');        		disp(Runtime');
disp('sigma: '); 			disp(sigmaP);
disp('b:');      			disp(bP);
disp('Tolerance Flag: ');		disp(ToleranceFlag);

figure(13);
title('C-1');
hold on;
plot(Timevec,fexact,'-ko','LineWidth',2);
plot(Timevec,real(fapprox),'-rx');
hold off;
xlabel('Time [sec]');
legend('Exact','Approximate');
drawnow;

saveas(gcf,[pwd,filesep,'EXAMPLE_FIGURES',filesep,'C_1.jpg']);

clearvars -except SetTest;
end

%%%%%%%%%%%%%%%
if(SetTest.C2==1)
% C) A-S_29.3.079
disp('------------------------------------------------------------');
disp('Test C-2');
disp('Wrapped Weeks Method: Alpha-Rho Search');

Fs = 'exp(1./s)./(s.^(1.5))';
ft = 'sinh(2*sqrt(t))/sqrt(pi)'
Timevec = 1:2:11;
RelErrorTol = [0.01 0.01 0.01 0.01 0.01 0.01]; %1 percent relative error

fexact = feval(inline(ft,'t'),Timevec);

tic;
[fapprox,alphaP,rhoP,EstRelError,EstAbsError,~,~,ToleranceFlag] = WeeksMethod(Fs,Timevec,RelErrorTol,2);
Runtime = toc;

MeaAbsError = 100*abs(fexact-fapprox);    %percent
MeaRelError = MeaAbsError./abs(fexact); %percent

disp('Time: ');  			disp(Timevec);
disp('Exact: '); 			disp(fexact); 
disp('Approximate: ');      		disp(real(fapprox)); 
disp('Measured Absolute Error %: ');  	disp(MeaAbsError); 
disp('Estimated Absolute Error %: '); 	disp(EstAbsError*100);
disp('Measured Relative Error %: ');  	disp(MeaRelError); 
disp('Estimated Relative Error %: '); 	disp(EstRelError*100);
disp('Run Time :');        		disp(Runtime');
disp('alpha: '); 			disp(alphaP);
disp('rho:');      			disp(rhoP);
disp('Tolerance Flag: ');		disp(ToleranceFlag);

figure(14);
title('C-2');
hold on;
plot(Timevec,fexact,'-ko','LineWidth',2);
plot(Timevec,real(fapprox),'-rx');
hold off;
xlabel('Time [sec]');
legend('Exact','Approximate');
drawnow;

saveas(gcf,[pwd,filesep,'EXAMPLE_FIGURES',filesep,'C_2.jpg']);

clearvars -except SetTest;
end

%%%%%%%%%%%%%%%
if(SetTest.C3==1)
% C) A-S_29.3.079
disp('------------------------------------------------------------');
disp('Test C-3');
disp('Wrapped Weeks Method: Adaptive Integrative alpha Fminsearch');

Fs = 'exp(1./s)./(s.^(1.5))';
ft = 'sinh(2*sqrt(t))/sqrt(pi)'
Timevec = 1:2:11;
RelErrorTol = [0.1 0.1 0.1 0.1 0.1 0.1]; %10 percent relative error

fexact = feval(inline(ft,'t'),Timevec);

tic;
[fapprox,alphaP,rhoP,EstRelError,EstAbsError,~,~,ToleranceFlag] = WeeksMethod(Fs,Timevec,RelErrorTol,5);
Runtime = toc;

MeaAbsError = 100*abs(fexact-fapprox);    %percent
MeaRelError = MeaAbsError./abs(fexact); %percent

disp('Time: ');  			disp(Timevec);
disp('Exact: '); 			disp(fexact); 
disp('Approximate: ');      		disp(real(fapprox)); 
disp('Measured Absolute Error %: ');  	disp(MeaAbsError); 
disp('Estimated Absolute Error %: '); 	disp(EstAbsError*100);
disp('Measured Relative Error %: ');  	disp(MeaRelError); 
disp('Estimated Relative Error %: '); 	disp(EstRelError*100);
disp('Run Time :');        		disp(Runtime');
disp('alpha: '); 			disp(alphaP);
disp('rho:');      			disp(rhoP);
disp('Tolerance Flag: ');		disp(ToleranceFlag);

figure(15);
title('C-3');
hold on;
plot(Timevec,fexact,'-ko','LineWidth',2);
plot(Timevec,real(fapprox),'-rx');
hold off;
xlabel('Time [sec]');
legend('Exact','Approximate');
drawnow;

saveas(gcf,[pwd,filesep,'EXAMPLE_FIGURES',filesep,'C_3.jpg']);

clearvars -except SetTest;
end


%%%%%%%%%%%%%%%
if(SetTest.D1==1)
% D) A-S_29.3.011 (challenging case for numerical inversion)
disp('------------------------------------------------------------');
disp('Test D-1');
disp('Wrapped Weeks Method: CPU Local Search');

Fs = '(s+0.6).^(-9/11)'
ft = '(t.^(-2/11))./(gamma(9/11)*exp(0.6*t))'
Timevec = 0.5:0.5:3;
RelErrorTol = [0.1 0.1 0.1 0.1 0.1 0.1]; %10 percent relative error for all times

fexact = feval(inline(ft,'t'),Timevec);

tic;
[fapprox,sigmaP,bP,EstRelError,EstAbsError,~,~,ToleranceFlag] = WeeksMethod(Fs,Timevec,RelErrorTol);
Runtime = toc;

MeaAbsError = 100*abs(fexact-fapprox);    %percent
MeaRelError = MeaAbsError./abs(fexact); %percent

disp('Time: ');  			disp(Timevec);
disp('Exact: '); 			disp(fexact); 
disp('Approximate: ');      		disp(real(fapprox)); 
disp('Measured Absolute Error %: ');  	disp(MeaAbsError); 
disp('Estimated Absolute Error %: '); 	disp(EstAbsError*100);
disp('Measured Relative Error %: ');  	disp(MeaRelError); 
disp('Estimated Relative Error %: '); 	disp(EstRelError*100);
disp('Run Time :');        		disp(Runtime');
disp('sigma: '); 			disp(sigmaP);
disp('b:');      			disp(bP);
disp('Tolerance Flag: ');		disp(ToleranceFlag);

figure(16);
title('D-1');
hold on;
plot(Timevec,fexact,'-ko','LineWidth',2);
plot(Timevec,real(fapprox),'-rx');
hold off;
xlabel('Time [sec]');
legend('Exact','Approximate');
drawnow;

saveas(gcf,[pwd,filesep,'EXAMPLE_FIGURES',filesep,'D_1.jpg']);

clearvars -except SetTest;
end

%%%%%%%%%%%%%%%
if(SetTest.D2==1)
% D) A-S_29.3.011 
disp('------------------------------------------------------------');
disp('Test D-2');
disp('Wrapped Weeks Method: Alpha-Rho Global Search');

Fs = '(s+0.6).^(-9/11)'
ft = '(t.^(-2/11))./(gamma(9/11)*exp(0.6*t))'
Timevec = 0.5:0.5:3;
RelErrorTol = [0.1 0.1 0.1 0.1 0.1 0.1]; %10 percent relative error for all times

fexact = feval(inline(ft,'t'),Timevec);

tic;
[fapprox,alphaP,rhoP,EstRelError,EstAbsError,~,~,ToleranceFlag] = WeeksMethod(Fs,Timevec,RelErrorTol,2);
Runtime = toc;

MeaAbsError = 100*abs(fexact-fapprox);    %percent
MeaRelError = MeaAbsError./abs(fexact); %percent

disp('Time: ');  			disp(Timevec);
disp('Exact: '); 			disp(fexact); 
disp('Approximate: ');      		disp(real(fapprox)); 
disp('Measured Absolute Error %: ');  	disp(MeaAbsError); 
disp('Estimated Absolute Error %: '); 	disp(EstAbsError*100);
disp('Measured Relative Error %: ');  	disp(MeaRelError); 
disp('Estimated Relative Error %: '); 	disp(EstRelError*100);
disp('Run Time :');        		disp(Runtime');
disp('alpha: '); 			disp(alphaP);
disp('rho:');      			disp(rhoP);
disp('Tolerance Flag: ');		disp(ToleranceFlag);

figure(17);
title('D-2');
hold on;
plot(Timevec,fexact,'-ko','LineWidth',2);
plot(Timevec,real(fapprox),'-rx');
hold off;
xlabel('Time [sec]');
legend('Exact','Approximate');
drawnow;

saveas(gcf,[pwd,filesep,'EXAMPLE_FIGURES',filesep,'D_2.jpg']);

clearvars -except SetTest;
end

%%%%%%%%%%%%%%%
if(SetTest.D3==1)
% D) A-S_29.3.011 
disp('------------------------------------------------------------');
disp('Test D-3');
disp('Wrapped Weeks Method: Adaptive Integrative alpha Fminsearch');

Fs = '(s+0.6).^(-9/11)'
ft = '(t.^(-2/11))./(gamma(9/11)*exp(0.6*t))'
Timevec = 0.5:0.5:3;
RelErrorTol = [0.1 0.1 0.1 0.1 0.1 0.1]; %10 percent relative error for all times

fexact = feval(inline(ft,'t'),Timevec);

tic;
[fapprox,alphaP,rhoP,EstRelError,EstAbsError,~,~,ToleranceFlag] = WeeksMethod(Fs,Timevec,RelErrorTol,5);
Runtime = toc;

MeaAbsError = 100*abs(fexact-fapprox);    %percent
MeaRelError = MeaAbsError./abs(fexact); %percent

disp('Time: ');  			disp(Timevec);
disp('Exact: '); 			disp(fexact); 
disp('Approximate: ');      		disp(real(fapprox)); 
disp('Measured Absolute Error %: ');  	disp(MeaAbsError); 
disp('Estimated Absolute Error %: '); 	disp(EstAbsError*100);
disp('Measured Relative Error %: ');  	disp(MeaRelError); 
disp('Estimated Relative Error %: '); 	disp(EstRelError*100);
disp('Run Time :');        		disp(Runtime');
disp('alpha: '); 			disp(alphaP);
disp('rho:');      			disp(rhoP);
disp('Tolerance Flag: ');		disp(ToleranceFlag);

figure(18);
title('D-3');
hold on;
plot(Timevec,fexact,'-ko','LineWidth',2);
plot(Timevec,real(fapprox),'-rx');
hold off;
xlabel('Time [sec]');
legend('Exact','Approximate');
drawnow;

saveas(gcf,[pwd,filesep,'EXAMPLE_FIGURES',filesep,'D_3.jpg']);

clearvars -except SetTest;
end

%%%%%%%%%%%%%%%
if(SetTest.E1==1)
 
disp('------------------------------------------------------------');
disp('Test E-1');
disp('Wrapped Weeks Method: Alpha-Rho fminsearch Search');

Fs = 'fnDuffyExample(s,3,3)' %(s,rho,r)
ft = '1'
Timevec = 0.5:0.5:3;
RelErrorTol = [0.1 0.1 0.1 0.1 0.1 0.1]; %10 percent relative error for all times

fexact = feval(inline(ft,'t'),Timevec);

tic;
[fapprox,alphaP,rhoP,EstRelError,EstAbsError,~,~,ToleranceFlag] = WeeksMethod(Fs,Timevec,RelErrorTol,2);
Runtime = toc;

MeaAbsError = 100*abs(fexact-fapprox);  %percent
MeaRelError = MeaAbsError./abs(fexact); %percent

disp('Time: ');                 disp(Timevec);
disp('Exact: ');                disp(fexact); 
disp('Approximate: ');      	disp(real(fapprox)); 
disp('Measured Absolute Error %: ');  	disp(MeaAbsError); 
disp('Estimated Absolute Error %: '); 	disp(EstAbsError*100);
disp('Measured Relative Error %: ');  	disp(MeaRelError); 
disp('Estimated Relative Error %: '); 	disp(EstRelError*100);
disp('Run Time :');        		disp(Runtime');
disp('alpha: ');                disp(alphaP);
disp('rho:');                   disp(rhoP);
disp('Tolerance Flag: ');		disp(ToleranceFlag);

figure(19);
title('E-1');
hold on;
plot(Timevec,fexact,'-ko','LineWidth',2);
plot(Timevec,real(fapprox),'-rx');
hold off;
xlabel('Time [sec]');
legend('Exact','Approximate');
drawnow;

saveas(gcf,[pwd,filesep,'EXAMPLE_FIGURES',filesep,'E_1.jpg']);

clearvars -except SetTest;
end

%%%%%%%%%%%%%%%
if(SetTest.E2==1)
 
disp('------------------------------------------------------------');
disp('Test E-2');
disp('Wrapped Weeks Method: Alpha-Rho fminbnd Search');

Fs = 'fnDuffyExample(s,3,3)' %(s,rho,r)
ft = '1'
Timevec = 0.5:0.5:3;
RelErrorTol = [0.001 0.001 0.001 0.001 0.001 0.001]; %0.1 percent relative error for all times

fexact = feval(inline(ft,'t'),Timevec);

tic;
[fapprox,alphaP,rhoP,EstRelError,EstAbsError,~,~,ToleranceFlag] = WeeksMethod(Fs,Timevec,RelErrorTol,3);
Runtime = toc;

MeaAbsError = 100*abs(fexact-fapprox);  %percent
MeaRelError = MeaAbsError./abs(fexact); %percent

disp('Time: ');                 disp(Timevec);
disp('Exact: ');                disp(fexact); 
disp('Approximate: ');      	disp(real(fapprox)); 
disp('Measured Absolute Error %: ');  	disp(MeaAbsError); 
disp('Estimated Absolute Error %: '); 	disp(EstAbsError*100);
disp('Measured Relative Error %: ');  	disp(MeaRelError); 
disp('Estimated Relative Error %: '); 	disp(EstRelError*100);
disp('Run Time :');        		disp(Runtime');
disp('alpha: ');                disp(alphaP);
disp('rho:');                   disp(rhoP);
disp('Tolerance Flag: ');		disp(ToleranceFlag);

figure(20);
title('E-2');
hold on;
plot(Timevec,fexact,'-ko','LineWidth',2);
plot(Timevec,real(fapprox),'-rx');
hold off;
xlabel('Time [sec]');
legend('Exact','Approximate');
drawnow;

saveas(gcf,[pwd,filesep,'EXAMPLE_FIGURES',filesep,'E_2.jpg']);

clearvars -except SetTest;
end

%%%%%%%%%%%%%%%
if(SetTest.E3==1)
disp('------------------------------------------------------------');
disp('Test E-3');
disp('Wrapped Weeks Method: Alpha-Rho Global Search');

Fs = 'fnDuffyExample(s,3,3)' %(s,rho,r)
ft = '1'
Timevec = 0.5:0.5:3;
RelErrorTol = [0.1 0.1 0.1 0.1 0.1 0.1]; %10 percent relative error for all times

fexact = feval(inline(ft,'t'),Timevec);

tic;
[fapprox,alphaP,rhoP,EstRelError,EstAbsError,~,~,ToleranceFlag] = WeeksMethod(Fs,Timevec,RelErrorTol,4);
Runtime = toc;

MeaAbsError = 100*abs(fexact-fapprox);  %percent
MeaRelError = MeaAbsError./abs(fexact); %percent

disp('Time: ');                 disp(Timevec);
disp('Exact: ');                disp(fexact); 
disp('Approximate: ');      	disp(real(fapprox)); 
disp('Measured Absolute Error %: ');  	disp(MeaAbsError); 
disp('Estimated Absolute Error %: '); 	disp(EstAbsError*100);
disp('Measured Relative Error %: ');  	disp(MeaRelError); 
disp('Estimated Relative Error %: '); 	disp(EstRelError*100);
disp('Run Time :');        		disp(Runtime');
disp('alpha: ');                disp(alphaP);
disp('rho:');                   disp(rhoP);
disp('Tolerance Flag: ');		disp(ToleranceFlag);

figure(21);
title('E-3');
hold on;
plot(Timevec,fexact,'-ko','LineWidth',2);
plot(Timevec,real(fapprox),'-rx');
hold off;
xlabel('Time [sec]');
legend('Exact','Approximate');
drawnow;

saveas(gcf,[pwd,filesep,'EXAMPLE_FIGURES',filesep,'E_3.jpg']);

clearvars -except SetTest;
end

%%%%%%%%%%%%%%%
if(SetTest.E4==1)
disp('------------------------------------------------------------');
disp('Test E-4');
disp('Wrapped Weeks Method: Alpha-Rho Global Search');

Fs = 'fnDuffyExample(s,3,3)' %(s,rho,r)
ft = '1'
Timevec = 0.5:0.5:3;
RelErrorTol = [0.1 0.1 0.1 0.1 0.1 0.1]; %10 percent relative error for all times

fexact = feval(inline(ft,'t'),Timevec);

tic;
[fapprox,alphaP,rhoP,EstRelError,EstAbsError,~,~,ToleranceFlag] = WeeksMethod(Fs,Timevec,RelErrorTol,5);
Runtime = toc;

MeaAbsError = 100*abs(fexact-fapprox);  %percent
MeaRelError = MeaAbsError./abs(fexact); %percent

disp('Time: ');                 disp(Timevec);
disp('Exact: ');                disp(fexact); 
disp('Approximate: ');      	disp(real(fapprox)); 
disp('Measured Absolute Error %: ');  	disp(MeaAbsError); 
disp('Estimated Absolute Error %: '); 	disp(EstAbsError*100);
disp('Measured Relative Error %: ');  	disp(MeaRelError); 
disp('Estimated Relative Error %: '); 	disp(EstRelError*100);
disp('Run Time :');        		disp(Runtime');
disp('alpha: ');                disp(alphaP);
disp('rho:');                   disp(rhoP);
disp('Tolerance Flag: ');		disp(ToleranceFlag);

figure(22);
title('E-4');
hold on;
plot(Timevec,fexact,'-ko','LineWidth',2);
plot(Timevec,real(fapprox),'-rx');
hold off;
xlabel('Time [sec]');
legend('Exact','Approximate');
drawnow;

saveas(gcf,[pwd,filesep,'EXAMPLE_FIGURES',filesep,'E_4.jpg']);

clearvars -except SetTest;
end


%End of Script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
