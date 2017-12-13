

set(0,'DefaultFigureVisible','off')
B = 75;
ep = 2;
cof = 0.5;
jmax = size(cof,2);
brommean = [];
weeksmean = [];


for j=1:jmax
tic()
T=51.3:10:301.3;
num = max(size(T));
sigma = cof(j)./T.^ep;
sigmaW = cof(j)./T;
b = B;

BromError = [];
WeeksError = [];


BromEstimate = zeros(1,num);
WeeksEstimate = zeros(1,num);

Bromtime = zeros(1,num);
Weekstime = zeros(1,num);

% D&M fcn 1

True = dm1(T);

for i=1:num
    fun = @(y)dml1(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-b,b,'AbsTol',1e-16);
    fun = 'dml1(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,sigmaW(i),0.5);
end

error = True-BromEstimate;
RMSE = sqrt(error*error'/num);
ABS = mean(abs(error));
REL = mean(abs(error./True));
BromError = [BromError;RMSE,ABS,REL];

error = True-WeeksEstimate;
RMSE = sqrt(error*error'/num);
ABS = mean(abs(error));
REL = mean(abs(error./True));
WeeksError = [WeeksError;RMSE,ABS,REL];


figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'o','MarkerSize',7,'Color','red')
plot(T,real(WeeksEstimate),'*','MarkerSize',7,'Color','blue')
title('DM 1')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
hold off

% D&M fcn 2

True = dm2(T);

for i=1:num
    fun = @(y)dml2(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-b,b,'AbsTol',1e-16);
    fun = 'dml2(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,sigmaW(i),0.5);
end

error = True-BromEstimate;
RMSE = sqrt(error*error'/num);
ABS = mean(abs(error));
REL = mean(abs(error./True));
BromError = [BromError;RMSE,ABS,REL];

error = True-WeeksEstimate;
RMSE = sqrt(error*error'/num);
ABS = mean(abs(error));
REL = mean(abs(error./True));
WeeksError = [WeeksError;RMSE,ABS,REL];

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'o','MarkerSize',7,'Color','red')
plot(T,real(WeeksEstimate),'*','MarkerSize',7,'Color','blue')
title('DM 2')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
hold off

% D&M fcn 5

True = dm5(T);

for i=1:num
    fun = @(y)dml5(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-b,b,'AbsTol',1e-16);
    fun = 'dml5(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,sigmaW(i),0.5);
end

error = True-BromEstimate;
RMSE = sqrt(error*error'/num);
ABS = mean(abs(error));
REL = mean(abs(error./True));
BromError = [BromError;RMSE,ABS,REL];

error = True-WeeksEstimate;
RMSE = sqrt(error*error'/num);
ABS = mean(abs(error));
REL = mean(abs(error./True));
WeeksError = [WeeksError;RMSE,ABS,REL];

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'o','MarkerSize',7,'Color','red')
plot(T,real(WeeksEstimate),'*','MarkerSize',7,'Color','blue')
title('DM 5')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
hold off

% D&M fcn 6

True = dm6(T);

for i=1:num
    fun = @(y)dml6(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-b,b,'AbsTol',1e-16);
    fun = 'dml6(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,sigmaW(i),0.5);
end

error = True-BromEstimate;
RMSE = sqrt(error*error'/num);
ABS = mean(abs(error));
REL = mean(abs(error./True));
BromError = [BromError;RMSE,ABS,REL];

error = True-WeeksEstimate;
RMSE = sqrt(error*error'/num);
ABS = mean(abs(error));
REL = mean(abs(error./True));
WeeksError = [WeeksError;RMSE,ABS,REL];

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'o','MarkerSize',7,'Color','red')
plot(T,real(WeeksEstimate),'*','MarkerSize',7,'Color','blue')
title('DM 6')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
hold off

% D&M fcn 8

True = dm8(T);

for i=1:num
    fun = @(y)dml8(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-b,b,'AbsTol',1e-16);
    fun = 'dml8(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,sigmaW(i),0.5);
end

error = True-BromEstimate;
RMSE = sqrt(error*error'/num);
ABS = mean(abs(error));
REL = mean(abs(error./True));
BromError = [BromError;RMSE,ABS,REL];

error = True-WeeksEstimate;
RMSE = sqrt(error*error'/num);
ABS = mean(abs(error));
REL = mean(abs(error./True));
WeeksError = [WeeksError;RMSE,ABS,REL];


figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'o','MarkerSize',7,'Color','red')
plot(T,real(WeeksEstimate),'*','MarkerSize',7,'Color','blue')
title('DM 8')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
hold off

% D&M fcn 9

True = dm9(T);

for i=1:num
    fun = @(y)dml9(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-b,b,'AbsTol',1e-16);
    fun = 'dml9(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,sigmaW(i),0.5);
end

error = True-BromEstimate;
RMSE = sqrt(error*error'/num);
ABS = mean(abs(error));
REL = mean(abs(error./True));
BromError = [BromError;RMSE,ABS,REL];

error = True-WeeksEstimate;
RMSE = sqrt(error*error'/num);
ABS = mean(abs(error));
REL = mean(abs(error./True));
WeeksError = [WeeksError;RMSE,ABS,REL];

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'o','MarkerSize',7,'Color','red')
plot(T,real(WeeksEstimate),'*','MarkerSize',7,'Color','blue')
title('DM 9')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
hold off

% D&M fcn 10

True = dm10(T);

for i=1:num
    fun = @(y)dml10(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-b,b,'AbsTol',1e-16);
    fun = 'dml10(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,sigmaW(i),0.5);
end

error = True-BromEstimate;
RMSE = sqrt(error*error'/num);
ABS = mean(abs(error));
REL = mean(abs(error./True));
BromError = [BromError;RMSE,ABS,REL];

error = True-WeeksEstimate;
RMSE = sqrt(error*error'/num);
ABS = mean(abs(error));
REL = mean(abs(error./True));
WeeksError = [WeeksError;RMSE,ABS,REL];

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'o','MarkerSize',7,'Color','red')
plot(T,real(WeeksEstimate),'*','MarkerSize',7,'Color','blue')
title('DM 10')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
hold off

% D&M fcn 11

True = dm11(T);

for i=1:num
    fun = @(y)dml11(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-b,b,'AbsTol',1e-16);
    fun = 'dml11(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,sigmaW(i),0.5);
end

error = True-BromEstimate;
RMSE = sqrt(error*error'/num);
ABS = mean(abs(error));
REL = mean(abs(error./True));
BromError = [BromError;RMSE,ABS,REL];

error = True-WeeksEstimate;
RMSE = sqrt(error*error'/num);
ABS = mean(abs(error));
REL = mean(abs(error./True));
WeeksError = [WeeksError;RMSE,ABS,REL];

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'o','MarkerSize',7,'Color','red')
plot(T,real(WeeksEstimate),'*','MarkerSize',7,'Color','blue')
title('DM 11')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
hold off

% % D&M fcn 12
% 
% True = dm12(T);
% 
% for i=1:num
%     fun = @(y)dml13(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
%     BromEstimate(i) = integral(fun,-b,b,'AbsTol',1e-16);
%     fun = 'dml12(s)';
%     WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,sigmaW(i),0.5);
% end
% 
% error = True-BromEstimate;
% RMSE = sqrt(error*error'/num);
% ABS = mean(abs(error));
% REL = mean(abs(error./True));
% BromError = [BromError;RMSE,ABS,REL];
% 
% error = True-WeeksEstimate;
% RMSE = sqrt(error*error'/num);
% ABS = mean(abs(error));
% REL = mean(abs(error./True));
% WeeksError = [WeeksError;RMSE,ABS,REL];
% 
% figure
% plot(T,True,'LineWidth',3,'Color','black')
% hold on
% plot(T,real(BromEstimate),'o','MarkerSize',7,'Color','red')
% plot(T,real(WeeksEstimate),'*','MarkerSize',7,'Color','blue')
% title('DM 12')
% xlabel('time')
% ylabel('f(t)')
% legend('True f(t)','Bromwich adapt','Weeks')
% hold off

% D&M fcn 13

True = dm13(T);

for i=1:num
    fun = @(y)dml13(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-b,b,'AbsTol',1e-16);
    fun = 'dml13(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,sigmaW(i),0.5);
end

error = True-BromEstimate;
RMSE = sqrt(error*error'/num);
ABS = mean(abs(error));
REL = mean(abs(error./True));
BromError = [BromError;RMSE,ABS,REL];

error = True-WeeksEstimate;
RMSE = sqrt(error*error'/num);
ABS = mean(abs(error));
REL = mean(abs(error./True));
WeeksError = [WeeksError;RMSE,ABS,REL];


figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'o','MarkerSize',7,'Color','red')
plot(T,real(WeeksEstimate),'*','MarkerSize',7,'Color','blue')
title('DM 13')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
hold off

% D&M fcn 15

True = dm15(T);

for i=1:num
    fun = @(y)dml15(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-b,b,'AbsTol',1e-16);
    fun = 'dml15(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,sigmaW(i),0.5);
end

error = True-BromEstimate;
RMSE = sqrt(error*error'/num);
ABS = mean(abs(error));
REL = mean(abs(error./True));
BromError = [BromError;RMSE,ABS,REL];

error = True-WeeksEstimate;
RMSE = sqrt(error*error'/num);
ABS = mean(abs(error));
REL = mean(abs(error./True));
WeeksError = [WeeksError;RMSE,ABS,REL];


figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'o','MarkerSize',7,'Color','red')
plot(T,real(WeeksEstimate),'*','MarkerSize',7,'Color','blue')
title('DM 15')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
hold off

% D&M fcn 16

True = dm16(T);

for i=1:num
    fun = @(y)dml16(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-b,b,'AbsTol',1e-16);
    fun = 'dml16(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,sigmaW(i),0.5);
end

error = True-BromEstimate;
RMSE = sqrt(error*error'/num);
ABS = mean(abs(error));
REL = mean(abs(error./True));
BromError = [BromError;RMSE,ABS,REL];

error = True-WeeksEstimate;
RMSE = sqrt(error*error'/num);
ABS = mean(abs(error));
REL = mean(abs(error./True));
WeeksError = [WeeksError;RMSE,ABS,REL];

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'o','MarkerSize',7,'Color','red')
plot(T,real(WeeksEstimate),'*','MarkerSize',7,'Color','blue')
title('DM 16')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
hold off

%A&V fcn 5

True = av5(T);

for i=1:num
    fun = @(y)avl5(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-b,b,'AbsTol',1e-16);
    fun = 'avl5(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,sigmaW(i),0.5);
end

error = True-BromEstimate;
RMSE = sqrt(error*error'/num);
ABS = mean(abs(error));
REL = mean(abs(error./True));
BromError = [BromError;RMSE,ABS,REL];

error = True-WeeksEstimate;
RMSE = sqrt(error*error'/num);
ABS = mean(abs(error));
REL = mean(abs(error./True));
WeeksError = [WeeksError;RMSE,ABS,REL];

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'o','MarkerSize',7,'Color','red')
plot(T,real(WeeksEstimate),'*','MarkerSize',7,'Color','blue')
title('AV 5')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
hold off

% A&V fcn 10

True = av10(T);

for i=1:num
    fun = @(y)avl10(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-b,b,'AbsTol',1e-16);
    fun = 'avl10(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,sigmaW(i),0.5);
end

error = True-BromEstimate;
RMSE = sqrt(error*error'/num);
ABS = mean(abs(error));
REL = mean(abs(error./True));
BromError = [BromError;RMSE,ABS,REL];

error = True-WeeksEstimate;
RMSE = sqrt(error*error'/num);
ABS = mean(abs(error));
REL = mean(abs(error./True));
WeeksError = [WeeksError;RMSE,ABS,REL];

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'o','MarkerSize',7,'Color','red')
plot(T,real(WeeksEstimate),'*','MarkerSize',7,'Color','blue')
title('AV 10')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
hold off



format long

weeksmean = mean(WeeksError)
brommean = [brommean;mean(BromError)]

toc()
end

brommean
weeksmean
