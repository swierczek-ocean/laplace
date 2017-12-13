T=201.2:1:300.2;
abcissa = 0.8./T;
abcissaw = 2./T;
b = 50;

BromError = [];
WeeksError = [];

BromEstimate = zeros(1,100);
WeeksEstimate = zeros(1,100);

Bromtime = zeros(1,100);
Weekstime = zeros(1,100);

% D&M fcn 1

True = dm1(T);

for i=1:100
    fun = @(x)dml1(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*b,abcissa(i)+1i*b,'AbsTol',1e-16);
    Bromtime(i) = toc();
    fun = 'dml1(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissaw(i),0.5);
    Weekstime(i) = toc();
end

error = True-BromEstimate;
RMSE = sqrt(error*error'/100);
ABS = mean(abs(error));
REL = mean(abs(error./True));
BromError = [BromError;RMSE,ABS,REL];

error = True-WeeksEstimate;
RMSE = sqrt(error*error'/100);
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

for i=1:100
    fun = @(x)dml2(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*b,abcissa(i)+1i*b,'AbsTol',1e-16);
    Bromtime(i) = toc();
    fun = 'dml2(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissaw(i),0.5);
    Weekstime(i) = toc();
end

error = True-BromEstimate;
RMSE = sqrt(error*error'/100);
ABS = mean(abs(error));
REL = mean(abs(error./True));
BromError = [BromError;RMSE,ABS,REL];

error = True-WeeksEstimate;
RMSE = sqrt(error*error'/100);
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

for i=1:100
    fun = @(x)dml5(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*b,abcissa(i)+1i*b,'AbsTol',1e-16);
    Bromtime(i) = toc();
    fun = 'dml5(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissaw(i),0.5);
    Weekstime(i) = toc();
end

error = True-BromEstimate;
RMSE = sqrt(error*error'/100);
ABS = mean(abs(error));
REL = mean(abs(error./True));
BromError = [BromError;RMSE,ABS,REL];

error = True-WeeksEstimate;
RMSE = sqrt(error*error'/100);
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

for i=1:100
    fun = @(x)dml6(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*b,abcissa(i)+1i*b,'AbsTol',1e-16);
    Bromtime(i) = toc();
    fun = 'dml6(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissaw(i),0.5);
    Weekstime(i) = toc();
end

error = True-BromEstimate;
RMSE = sqrt(error*error'/100);
ABS = mean(abs(error));
REL = mean(abs(error./True));
BromError = [BromError;RMSE,ABS,REL];

error = True-WeeksEstimate;
RMSE = sqrt(error*error'/100);
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

for i=1:100
    fun = @(x)dml8(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*b,abcissa(i)+1i*b);
    Bromtime(i) = toc();
    fun = 'dml8(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissaw(i),0.5);
    Weekstime(i) = toc();
end

error = True-BromEstimate;
RMSE = sqrt(error*error'/100);
ABS = mean(abs(error));
REL = mean(abs(error./True));
BromError = [BromError;RMSE,ABS,REL];

error = True-WeeksEstimate;
RMSE = sqrt(error*error'/100);
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

for i=1:100
    fun = @(x)dml9(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*b,abcissa(i)+1i*b,'AbsTol',1e-16);
    Bromtime(i) = toc();
    fun = 'dml9(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissaw(i),0.5);
    Weekstime(i) = toc();
end

error = True-BromEstimate;
RMSE = sqrt(error*error'/100);
ABS = mean(abs(error));
REL = mean(abs(error./True));
BromError = [BromError;RMSE,ABS,REL];

error = True-WeeksEstimate;
RMSE = sqrt(error*error'/100);
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

for i=1:100
    fun = @(x)dml10(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*b,abcissa(i)+1i*b,'AbsTol',1e-16);
    Bromtime(i) = toc();
    fun = 'dml10(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissaw(i),0.5);
    Weekstime(i) = toc();
end

error = True-BromEstimate;
RMSE = sqrt(error*error'/100);
ABS = mean(abs(error));
REL = mean(abs(error./True));
BromError = [BromError;RMSE,ABS,REL];

error = True-WeeksEstimate;
RMSE = sqrt(error*error'/100);
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

for i=1:100
    fun = @(x)dml11(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*b,abcissa(i)+1i*b,'AbsTol',1e-16);
    Bromtime(i) = toc();
    fun = 'dml11(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissaw(i),0.5);
    Weekstime(i) = toc();
end

error = True-BromEstimate;
RMSE = sqrt(error*error'/100);
ABS = mean(abs(error));
REL = mean(abs(error./True));
BromError = [BromError;RMSE,ABS,REL];

error = True-WeeksEstimate;
RMSE = sqrt(error*error'/100);
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

% D&M fcn 12

True = dm12(T);

for i=1:100
    fun = @(x)dml12(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*b,abcissa(i)+1i*b,'AbsTol',1e-16);
    Bromtime(i) = toc();
    fun = 'dml12(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissaw(i),0.5);
    Weekstime(i) = toc();
end

error = True-BromEstimate;
RMSE = sqrt(error*error'/100);
ABS = mean(abs(error));
REL = mean(abs(error./True));
BromError = [BromError;RMSE,ABS,REL];

error = True-WeeksEstimate;
RMSE = sqrt(error*error'/100);
ABS = mean(abs(error));
REL = mean(abs(error./True));
WeeksError = [WeeksError;RMSE,ABS,REL];

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'o','MarkerSize',7,'Color','red')
plot(T,real(WeeksEstimate),'*','MarkerSize',7,'Color','blue')
title('DM 12')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
hold off

% D&M fcn 13

True = dm13(T);

for i=1:100
    fun = @(x)dml13(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*b,abcissa(i)+1i*b,'AbsTol',1e-16);
    Bromtime(i) = toc();
    fun = 'dml13(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissaw(i),0.5);
    Weekstime(i) = toc();
end

error = True-BromEstimate;
RMSE = sqrt(error*error'/100);
ABS = mean(abs(error));
REL = mean(abs(error./True));
BromError = [BromError;RMSE,ABS,REL];

error = True-WeeksEstimate;
RMSE = sqrt(error*error'/100);
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

for i=1:100
    fun = @(x)dml15(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*b,abcissa(i)+1i*b,'AbsTol',1e-16);
    Bromtime(i) = toc();
    fun = 'dml15(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissaw(i),0.5);
    Weekstime(i) = toc();
end

error = True-BromEstimate;
RMSE = sqrt(error*error'/100);
ABS = mean(abs(error));
REL = mean(abs(error./True));
BromError = [BromError;RMSE,ABS,REL];

error = True-WeeksEstimate;
RMSE = sqrt(error*error'/100);
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

for i=1:100
    fun = @(x)dml16(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*b,abcissa(i)+1i*b,'AbsTol',1e-16);
    Bromtime(i) = toc();
    fun = 'dml16(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissaw(i),0.5);
    Weekstime(i) = toc();
end

error = True-BromEstimate;
RMSE = sqrt(error*error'/100);
ABS = mean(abs(error));
REL = mean(abs(error./True));
BromError = [BromError;RMSE,ABS,REL];

error = True-WeeksEstimate;
RMSE = sqrt(error*error'/100);
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

for i=1:100
    fun = @(x)avl5(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*b,abcissa(i)+1i*b,'AbsTol',1e-16);
    Bromtime(i) = toc();
    fun = 'avl5(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissaw(i),0.5);
    Weekstime(i) = toc();
end

error = True-BromEstimate;
RMSE = sqrt(error*error'/100);
ABS = mean(abs(error));
REL = mean(abs(error./True));
BromError = [BromError;RMSE,ABS,REL];

error = True-WeeksEstimate;
RMSE = sqrt(error*error'/100);
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

for i=1:100
    fun = @(x)avl10(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*b,abcissa(i)+1i*b,'AbsTol',1e-16);
    Bromtime(i) = toc();
    fun = 'avl10(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissaw(i),0.5);
    Weekstime(i) = toc();
end

error = True-BromEstimate;
RMSE = sqrt(error*error'/100);
ABS = mean(abs(error));
REL = mean(abs(error./True));
BromError = [BromError;RMSE,ABS,REL];

error = True-WeeksEstimate;
RMSE = sqrt(error*error'/100);
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

BromError
WeeksError





