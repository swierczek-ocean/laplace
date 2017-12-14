tic()
set(0,'DefaultFigureVisible','on')
B = 20;
ep = 2;
cof = 0.5;
jmax = size(cof,2);
brommean = [];
weeksmean = [];
a=1;

tic()
T=0.01;
num = max(size(T));
sigma = cof./T.^ep;
sigmaW = cof./T;

% sigma = ones(num);
% sigmaW = sigma;
b=0.5;

BromError = [];
WeeksError = [];


BromEstimate = zeros(1,num);
WeeksEstimate = zeros(1,num);

Bromtime = zeros(1,num);
Weekstime = zeros(1,num);

% S fcn 1

True = ilapt1(T);

for i=1:num
    fun = @(y)lapt1(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
    fun = 'lapt1(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,sigmaW(i),b);
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
title('ILT 1')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT1','-djpeg')
hold off


% S fcn 2

True = ilapt2(T);

for i=1:num
    fun = @(y)lapt2(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
    fun = 'lapt2(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,sigmaW(i),b);
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
title('ILT 2')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT2','-djpeg')
hold off

% S fcn 3

True = ilapt3(T);

for i=1:num
    fun = @(y)lapt3(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
    fun = 'lapt3(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,sigmaW(i),b);
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
title('ILT 3')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT3','-djpeg')
hold off

% S fcn 4

True = ilapt4(T);

for i=1:num
    fun = @(y)lapt4(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
    fun = 'lapt4(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,sigmaW(i),b);
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
title('ILT 4')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT4','-djpeg')
hold off

% S fcn 5

True = ilapt5(T);

for i=1:num
    fun = @(y)lapt5(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-11);
    fun = 'lapt5(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,sigmaW(i),b);
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
title('ILT 5')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT5','-djpeg')
hold off


% % S fcn 7

True = ilapt7(T);

for i=1:num
    fun = @(y)lapt7(a+ sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
    fun = 'lapt7(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,a+sigmaW(i),b);
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
title('ILT 7')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT7','-djpeg')
hold off


% S fcn 8

True = ilapt8(T);

for i=1:num
    fun = @(y)lapt8(a+sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
    fun = 'lapt8(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,a+ sigmaW(i),b);
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
title('ILT 8')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT8','-djpeg')
hold off

% S fcn 9

True = ilapt9(T);

for i=1:num
    fun = @(y)lapt9(a+sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
    fun = 'lapt9(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,a+sigmaW(i),b);
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
title('ILT 9')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT9','-djpeg')
hold off

% S fcn 10

True = ilapt10(T);

for i=1:num
    fun = @(y)lapt10(a+sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
    fun = 'lapt10(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,a+sigmaW(i),b);
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
title('ILT 10')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT10','-djpeg')
hold off

% S fcn 11

True = ilapt11(T);

for i=1:num
    fun = @(y)lapt11(a+sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-11);
    fun = 'lapt11(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,a+sigmaW(i),b);
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
title('ILT 11')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT11','-djpeg')
hold off


% S fcn 12

True = ilapt12(T);

for i=1:num
    fun = @(y)lapt12(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
    fun = 'lapt12(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,sigmaW(i),b);
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
title('ILT 12')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT12','-djpeg')
hold off


% S fcn 13

True = ilapt13(T);

for i=1:num
    fun = @(y)lapt13(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
    fun = 'lapt13(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,sigmaW(i),b);
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
title('ILT 13')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT13','-djpeg')
hold off

% % S fcn 14

True = ilapt14(T);

for i=1:num
    fun = @(y)lapt14(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
    fun = 'lapt14(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,sigmaW(i),b);
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
title('ILT 14')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT14','-djpeg')
hold off

% S fcn 15

True = ilapt15(T);

for i=1:num
    fun = @(y)lapt15(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
    fun = 'lapt15(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,sigmaW(i),b);
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
title('ILT 15')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT15','-djpeg')
hold off

% S fcn 16

True = ilapt16(T);

for i=1:num
    fun = @(y)lapt16(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-11);
    fun = 'lapt16(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,sigmaW(i),b);
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
title('ILT 16')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT16','-djpeg')
hold off


% S fcn 17

True = ilapt17(T);

for i=1:num
    fun = @(y)lapt17(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
    fun = 'lapt17(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,sigmaW(i),b);
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
title('ILT 17')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT17','-djpeg')
hold off


% S fcn 18

True = ilapt18(T);

for i=1:num
    fun = @(y)lapt18(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
    fun = 'lapt18(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,sigmaW(i),b);
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
title('ILT 18')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT18','-djpeg')
hold off

% S fcn 19

True = ilapt19(T);

for i=1:num
    fun = @(y)lapt19(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
    fun = 'lapt19(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,sigmaW(i),b);
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
title('ILT 19')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT19','-djpeg')
hold off

% S fcn 20

True = ilapt20(T);

for i=1:num
    fun = @(y)lapt20(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
    fun = 'lapt20(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,sigmaW(i),b);
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
title('ILT 20')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT20','-djpeg')
hold off

% % S fcn 21

True = ilapt21(T);

for i=1:num
    fun = @(y)lapt21(a+sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-11);
    fun = 'lapt21(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,a+sigmaW(i),b);
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
title('ILT 21')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT21','-djpeg')
hold off


% S fcn 22

True = ilapt22(T);

for i=1:num
    fun = @(y)lapt22(a+sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
    fun = 'lapt22(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,a+sigmaW(i),b);
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
title('ILT 22')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT22','-djpeg')
hold off


% S fcn 23

True = ilapt23(T);

for i=1:num
    fun = @(y)lapt23(a+sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
    fun = 'lapt23(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,a+sigmaW(i),b);
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
title('ILT 23')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT23','-djpeg')
hold off

% S fcn 24

True = ilapt24(T);

for i=1:num
    fun = @(y)lapt24(a+sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
    fun = 'lapt24(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,a+sigmaW(i),b);
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
title('ILT 24')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT24','-djpeg')
hold off

% S fcn 25

True = ilapt25(T);

for i=1:num
    fun = @(y)lapt25(a+sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
    fun = 'lapt25(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,a+sigmaW(i),b);
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
title('ILT 25')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT25','-djpeg')
hold off

% S fcn 26

True = ilapt26(T);

for i=1:num
    fun = @(y)lapt26(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-11);
    fun = 'lapt26(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,sigmaW(i),b);
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
title('ILT 26')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT26','-djpeg')
hold off


% S fcn 27

True = ilapt27(T);

for i=1:num
    fun = @(y)lapt27(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
    fun = 'lapt27(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,sigmaW(i),b);
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
title('ILT 27')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT27','-djpeg')
hold off


% S fcn 28

True = ilapt28(T);

for i=1:num
    fun = @(y)lapt28(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
    fun = 'lapt28(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,sigmaW(i),b);
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
title('ILT 28')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT28','-djpeg')
hold off

% S fcn 29

True = ilapt29(T);

for i=1:num
    fun = @(y)lapt29(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
    fun = 'lapt29(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,sigmaW(i),b);
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
title('ILT 29')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT29','-djpeg')
hold off

% S fcn 30

True = ilapt30(T);

for i=1:num
    fun = @(y)lapt30(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
    fun = 'lapt30(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,sigmaW(i),b);
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
title('ILT 30')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT30','-djpeg')
hold off

% % S fcn 31

True = ilapt31(T);

for i=1:num
    fun = @(y)lapt31(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-11);
    fun = 'lapt31(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,sigmaW(i),b);
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
title('ILT 31')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT31','-djpeg')
hold off


% S fcn 32

True = ilapt32(T);

for i=1:num
    fun = @(y)lapt32(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
    fun = 'lapt32(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,sigmaW(i),b);
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
title('ILT 32')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT32','-djpeg')
hold off


% S fcn 33

True = ilapt33(T);

for i=1:num
    fun = @(y)lapt33(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
    fun = 'lapt33(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,sigmaW(i),b);
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
title('ILT 33')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT33','-djpeg')
hold off

% S fcn 34

True = ilapt34(T);

for i=1:num
    fun = @(y)lapt34(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
    fun = 'lapt34(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,sigmaW(i),b);
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
title('ILT 34')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT34','-djpeg')
hold off

% S fcn 35

True = ilapt35(T);

for i=1:num
    fun = @(y)lapt35(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
    fun = 'lapt35(s)';
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),512,sigmaW(i),b);
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
title('ILT 35')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('ILT35','-djpeg')
hold off







format long

BromError
WeeksError

weeksmean = mean(WeeksError)
brommean = mean(BromError)

toc()






















