T=1:1:100;
abcissa = 1./T;

BromError = zeros(11,5);
WeeksError = zeros(11,5);

BromEstimate = zeros(1,100);
WeeksEstimate = zeros(1,100);

Bromtime = zeros(1,100);
Weekstime = zeros(1,100);

% D&M fcn 2

True = dm2(T);

for i=1:100
    fun = @(x)dml2(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*10,abcissa(i)+1i*10);
    Bromtime(i) = toc();
    fun = 'dml2(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissa(i),0.5);
    Weekstime(i) = toc();
end

error = True-BromEstimate;
BromError(1,:) = [sum(Bromtime),max(Bromtime),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-WeeksEstimate;
WeeksError(1,:) = [sum(Weekstime),max(Weekstime),mean(abs(error)),norm(error,2),max(abs(error))];


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
    BromEstimate(i) = integral(fun,abcissa(i)-1i*10,abcissa(i)+1i*10);
    Bromtime(i) = toc();
    fun = 'dml5(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissa(i),0.5);
    Weekstime(i) = toc();
end

error = True-BromEstimate;
BromError(2,:) = [sum(Bromtime),max(Bromtime),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-WeeksEstimate;
WeeksError(2,:) = [sum(Weekstime),max(Weekstime),mean(abs(error)),norm(error,2),max(abs(error))];

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

% D&M fcn 6

True = dm6(T);

for i=1:100
    fun = @(x)dml6(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*10,abcissa(i)+1i*10);
    Bromtime(i) = toc();
    fun = 'dml6(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissa(i),0.5);
    Weekstime(i) = toc();
end

error = True-BromEstimate;
BromError(3,:) = [sum(Bromtime),max(Bromtime),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-WeeksEstimate;
WeeksError(3,:) = [sum(Weekstime),max(Weekstime),mean(abs(error)),norm(error,2),max(abs(error))];

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

% D&M fcn 9

True = dm9(T);

for i=1:100
    fun = @(x)dml9(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*10,abcissa(i)+1i*10);
    Bromtime(i) = toc();
    fun = 'dml9(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissa(i),0.5);
    Weekstime(i) = toc();
end

error = True-BromEstimate;
BromError(4,:) = [sum(Bromtime),max(Bromtime),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-WeeksEstimate;
WeeksError(4,:) = [sum(Weekstime),max(Weekstime),mean(abs(error)),norm(error,2),max(abs(error))];

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
    BromEstimate(i) = integral(fun,abcissa(i)-1i*10,abcissa(i)+1i*10);
    Bromtime(i) = toc();
    fun = 'dml10(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissa(i),0.5);
    Weekstime(i) = toc();
end

error = True-BromEstimate;
BromError(5,:) = [sum(Bromtime),max(Bromtime),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-WeeksEstimate;
WeeksError(5,:) = [sum(Weekstime),max(Weekstime),mean(abs(error)),norm(error,2),max(abs(error))];

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
    BromEstimate(i) = integral(fun,abcissa(i)-1i*10,abcissa(i)+1i*10);
    Bromtime(i) = toc();
    fun = 'dml11(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissa(i),0.5);
    Weekstime(i) = toc();
end

error = True-BromEstimate;
BromError(6,:) = [sum(Bromtime),max(Bromtime),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-WeeksEstimate;
WeeksError(6,:) = [sum(Weekstime),max(Weekstime),mean(abs(error)),norm(error,2),max(abs(error))];

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
    BromEstimate(i) = integral(fun,abcissa(i)-1i*10,abcissa(i)+1i*10);
    Bromtime(i) = toc();
    fun = 'dml12(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissa(i),0.5);
    Weekstime(i) = toc();
end

error = True-BromEstimate;
BromError(7,:) = [sum(Bromtime),max(Bromtime),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-WeeksEstimate;
WeeksError(7,:) = [sum(Weekstime),max(Weekstime),mean(abs(error)),norm(error,2),max(abs(error))];

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

% D&M fcn 16

True = dm16(T);

for i=1:100
    fun = @(x)dml16(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*10,abcissa(i)+1i*10);
    Bromtime(i) = toc();
    fun = 'dml16(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissa(i),0.5);
    Weekstime(i) = toc();
end

error = True-BromEstimate;
BromError(8,:) = [sum(Bromtime),max(Bromtime),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-WeeksEstimate;
WeeksError(8,:) = [sum(Weekstime),max(Weekstime),mean(abs(error)),norm(error,2),max(abs(error))];

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

% A&V fcn 5

True = av5(T);

for i=1:100
    fun = @(x)avl5(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*10,abcissa(i)+1i*10);
    Bromtime(i) = toc();
    fun = 'avl5(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissa(i),0.5);
    Weekstime(i) = toc();
end

error = True-BromEstimate;
BromError(9,:) = [sum(Bromtime),max(Bromtime),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-WeeksEstimate;
WeeksError(9,:) = [sum(Weekstime),max(Weekstime),mean(abs(error)),norm(error,2),max(abs(error))];

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
    BromEstimate(i) = integral(fun,abcissa(i)-1i*10,abcissa(i)+1i*10);
    Bromtime(i) = toc();
    fun = 'avl10(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissa(i),0.5);
    Weekstime(i) = toc();
end

error = True-BromEstimate;
BromError(10,:) = [sum(Bromtime),max(Bromtime),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-WeeksEstimate;
WeeksError(10,:) = [sum(Weekstime),max(Weekstime),mean(abs(error)),norm(error,2),max(abs(error))];

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

% A&V fcn 12

True = av12(T);

for i=1:100
    fun = @(x)avl12(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*10,abcissa(i)+1i*10);
    Bromtime(i) = toc();
    fun = 'avl12(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissa(i),0.5);
    Weekstime(i) = toc();
end


error = True-BromEstimate;
BromError(11,:) = [sum(Bromtime),max(Bromtime),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-WeeksEstimate;
WeeksError(11,:) = [sum(Weekstime),max(Weekstime),mean(abs(error)),norm(error,2),max(abs(error))];

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'o','MarkerSize',7,'Color','red')
plot(T,real(WeeksEstimate),'*','MarkerSize',7,'Color','blue')
title('AV 12')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
hold off

format long g

BromError
WeeksError