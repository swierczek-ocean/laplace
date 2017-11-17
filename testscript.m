T=1:1:100;
abcissa = 1./T;

Wrap1Error = zeros(11,5);
Wrap2Error = zeros(11,5);
BromError = zeros(11,5);
WeeksError = zeros(11,5);

Wrap1Estimate = zeros(1,100);
Wrap2Estimate = zeros(1,100);
BromEstimate = zeros(1,100);
WeeksEstimate = zeros(1,100);

Wrap1time = zeros(1,100);
Wrap2time = zeros(1,100);
Bromtime = zeros(1,100);
Weekstime = zeros(1,100);

% D&M fcn 2

True = dm2(T);

for i=1:100
    fun = @(x)dml2(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*10,abcissa(i)+1i*10);
    Bromtime(i) = toc();
    tic();
    Wrap1Estimate(i) = contourint(1000,abcissa(i),fun);
    Wrap1time(i) = toc();
    tic();
    Wrap2Estimate(i) = contourint2(-0.1,0+abcissa(i),-0.1,0.1,.0001,fun);
    Wrap2time(i) = toc();
    fun = 'dml2(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissa(i),0.5);
    Weekstime(i) = toc();
end

error = True-Wrap1Estimate;
Wrap1Error(1,:) = [sum(Wrap1time),max(Wrap1time),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-Wrap2Estimate;
Wrap2Error(1,:) = [sum(Wrap2time),max(Wrap2time),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-BromEstimate;
BromError(1,:) = [sum(Bromtime),max(Bromtime),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-WeeksEstimate;
WeeksError(1,:) = [sum(Weekstime),max(Weekstime),mean(abs(error)),norm(error,2),max(abs(error))];


figure
plot(T,True,'LineWidth',4,'Color','black')
hold on
plot(T,real(Wrap1Estimate),'o','MarkerSize',7,'Color','green')
plot(T,real(Wrap2Estimate),'o','MarkerSize',7,'Color','cyan')
plot(T,real(BromEstimate),'o','MarkerSize',7,'Color','red')
plot(T,real(WeeksEstimate),'o','MarkerSize',7,'Color','blue')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','circular trapez','square trapez','Bromwich adapt','Weeks')
hold off

% D&M fcn 5

True = dm5(T);

for i=1:100
    fun = @(x)dml5(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*10,abcissa(i)+1i*10);
    Bromtime(i) = toc();
    tic();
    Wrap1Estimate(i) = contourint(1000,abcissa(i),fun);
    Wrap1time(i) = toc();
    tic();
    Wrap2Estimate(i) = contourint2(-0.1,0+abcissa(i),-0.1,0.1,.0001,fun);
    Wrap2time(i) = toc();
    fun = 'dml5(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissa(i),0.5);
    Weekstime(i) = toc();
end

error = True-Wrap1Estimate;
Wrap1Error(2,:) = [sum(Wrap1time),max(Wrap1time),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-Wrap2Estimate;
Wrap2Error(2,:) = [sum(Wrap2time),max(Wrap2time),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-BromEstimate;
BromError(2,:) = [sum(Bromtime),max(Bromtime),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-WeeksEstimate;
WeeksError(2,:) = [sum(Weekstime),max(Weekstime),mean(abs(error)),norm(error,2),max(abs(error))];


figure
plot(T,True,'LineWidth',4,'Color','black')
hold on
plot(T,real(Wrap1Estimate),'o','MarkerSize',7,'Color','green')
plot(T,real(Wrap2Estimate),'o','MarkerSize',7,'Color','cyan')
plot(T,real(BromEstimate),'o','MarkerSize',7,'Color','red')
plot(T,real(WeeksEstimate),'o','MarkerSize',7,'Color','blue')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','circular trapez','square trapez','Bromwich adapt','Weeks')
hold off

% D&M fcn 6

True = dm6(T);

for i=1:100
    fun = @(x)dml6(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*10,abcissa(i)+1i*10);
    Bromtime(i) = toc();
    tic();
    Wrap1Estimate(i) = contourint(1000,abcissa(i),fun);
    Wrap1time(i) = toc();
    tic();
    Wrap2Estimate(i) = contourint2(-0.1,0+abcissa(i),-0.1,0.1,.0001,fun);
    Wrap2time(i) = toc();
    fun = 'dml6(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissa(i),0.5);
    Weekstime(i) = toc();
end

error = True-Wrap1Estimate;
Wrap1Error(3,:) = [sum(Wrap1time),max(Wrap1time),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-Wrap2Estimate;
Wrap2Error(3,:) = [sum(Wrap2time),max(Wrap2time),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-BromEstimate;
BromError(3,:) = [sum(Bromtime),max(Bromtime),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-WeeksEstimate;
WeeksError(3,:) = [sum(Weekstime),max(Weekstime),mean(abs(error)),norm(error,2),max(abs(error))];

figure
plot(T,True,'LineWidth',4,'Color','black')
hold on
plot(T,real(Wrap1Estimate),'o','MarkerSize',7,'Color','green')
plot(T,real(Wrap2Estimate),'o','MarkerSize',7,'Color','cyan')
plot(T,real(BromEstimate),'o','MarkerSize',7,'Color','red')
plot(T,real(WeeksEstimate),'o','MarkerSize',7,'Color','blue')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','circular trapez','square trapez','Bromwich adapt','Weeks')
hold off

% D&M fcn 9

True = dm9(T);

for i=1:100
    fun = @(x)dml9(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*10,abcissa(i)+1i*10);
    Bromtime(i) = toc();
    tic();
    Wrap1Estimate(i) = contourint(1000,abcissa(i),fun);
    Wrap1time(i) = toc();
    tic();
    Wrap2Estimate(i) = contourint2(-0.1,0+abcissa(i),-0.1,0.1,.0001,fun);
    Wrap2time(i) = toc();
    fun = 'dml9(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissa(i),0.5);
    Weekstime(i) = toc();
end

error = True-Wrap1Estimate;
Wrap1Error(4,:) = [sum(Wrap1time),max(Wrap1time),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-Wrap2Estimate;
Wrap2Error(4,:) = [sum(Wrap2time),max(Wrap2time),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-BromEstimate;
BromError(4,:) = [sum(Bromtime),max(Bromtime),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-WeeksEstimate;
WeeksError(4,:) = [sum(Weekstime),max(Weekstime),mean(abs(error)),norm(error,2),max(abs(error))];

figure
plot(T,True,'LineWidth',4,'Color','black')
hold on
plot(T,real(Wrap1Estimate),'o','MarkerSize',7,'Color','green')
plot(T,real(Wrap2Estimate),'o','MarkerSize',7,'Color','cyan')
plot(T,real(BromEstimate),'o','MarkerSize',7,'Color','red')
plot(T,real(WeeksEstimate),'o','MarkerSize',7,'Color','blue')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','circular trapez','square trapez','Bromwich adapt','Weeks')
hold off

% D&M fcn 10

True = dm10(T);

for i=1:100
    fun = @(x)dml10(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*10,abcissa(i)+1i*10);
    Bromtime(i) = toc();
    tic();
    Wrap1Estimate(i) = contourint(1000,abcissa(i),fun);
    Wrap1time(i) = toc();
    tic();
    Wrap2Estimate(i) = contourint2(-0.1,0+abcissa(i),-0.1,0.1,.0001,fun);
    Wrap2time(i) = toc();
    fun = 'dml10(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissa(i),0.5);
    Weekstime(i) = toc();
end

error = True-Wrap1Estimate;
Wrap1Error(5,:) = [sum(Wrap1time),max(Wrap1time),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-Wrap2Estimate;
Wrap2Error(5,:) = [sum(Wrap2time),max(Wrap2time),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-BromEstimate;
BromError(5,:) = [sum(Bromtime),max(Bromtime),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-WeeksEstimate;
WeeksError(5,:) = [sum(Weekstime),max(Weekstime),mean(abs(error)),norm(error,2),max(abs(error))];

figure
plot(T,True,'LineWidth',4,'Color','black')
hold on
plot(T,real(Wrap1Estimate),'o','MarkerSize',7,'Color','green')
plot(T,real(Wrap2Estimate),'o','MarkerSize',7,'Color','cyan')
plot(T,real(BromEstimate),'o','MarkerSize',7,'Color','red')
plot(T,real(WeeksEstimate),'o','MarkerSize',7,'Color','blue')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','circular trapez','square trapez','Bromwich adapt','Weeks')
hold off

% D&M fcn 11

True = dm11(T);

for i=1:100
    fun = @(x)dml11(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*10,abcissa(i)+1i*10);
    Bromtime(i) = toc();
    tic();
    Wrap1Estimate(i) = contourint(1000,abcissa(i),fun);
    Wrap1time(i) = toc();
    tic();
    Wrap2Estimate(i) = contourint2(-0.1,0+abcissa(i),-0.1,0.1,.0001,fun);
    Wrap2time(i) = toc();
    fun = 'dml11(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissa(i),0.5);
    Weekstime(i) = toc();
end

error = True-Wrap1Estimate;
Wrap1Error(6,:) = [sum(Wrap1time),max(Wrap1time),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-Wrap2Estimate;
Wrap2Error(6,:) = [sum(Wrap2time),max(Wrap2time),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-BromEstimate;
BromError(6,:) = [sum(Bromtime),max(Bromtime),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-WeeksEstimate;
WeeksError(6,:) = [sum(Weekstime),max(Weekstime),mean(abs(error)),norm(error,2),max(abs(error))];

figure
plot(T,True,'LineWidth',4,'Color','black')
hold on
plot(T,real(Wrap1Estimate),'o','MarkerSize',7,'Color','green')
plot(T,real(Wrap2Estimate),'o','MarkerSize',7,'Color','cyan')
plot(T,real(BromEstimate),'o','MarkerSize',7,'Color','red')
plot(T,real(WeeksEstimate),'o','MarkerSize',7,'Color','blue')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','circular trapez','square trapez','Bromwich adapt','Weeks')
hold off

% D&M fcn 12

True = dm12(T);

for i=1:100
    fun = @(x)dml12(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*10,abcissa(i)+1i*10);
    Bromtime(i) = toc();
    tic();
    Wrap1Estimate(i) = contourint(1000,abcissa(i),fun);
    Wrap1time(i) = toc();
    tic();
    Wrap2Estimate(i) = contourint2(-0.1,0+abcissa(i),-0.1,0.1,.0001,fun);
    Wrap2time(i) = toc();
    fun = 'dml12(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissa(i),0.5);
    Weekstime(i) = toc();
end

error = True-Wrap1Estimate;
Wrap1Error(7,:) = [sum(Wrap1time),max(Wrap1time),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-Wrap2Estimate;
Wrap2Error(7,:) = [sum(Wrap2time),max(Wrap2time),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-BromEstimate;
BromError(7,:) = [sum(Bromtime),max(Bromtime),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-WeeksEstimate;
WeeksError(7,:) = [sum(Weekstime),max(Weekstime),mean(abs(error)),norm(error,2),max(abs(error))];

figure
plot(T,True,'LineWidth',4,'Color','black')
hold on
plot(T,real(Wrap1Estimate),'o','MarkerSize',7,'Color','green')
plot(T,real(Wrap2Estimate),'o','MarkerSize',7,'Color','cyan')
plot(T,real(BromEstimate),'o','MarkerSize',7,'Color','red')
plot(T,real(WeeksEstimate),'o','MarkerSize',7,'Color','blue')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','circular trapez','square trapez','Bromwich adapt','Weeks')
hold off

% D&M fcn 16

True = dm16(T);

for i=1:100
    fun = @(x)dml16(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*10,abcissa(i)+1i*10);
    Bromtime(i) = toc();
    tic();
    Wrap1Estimate(i) = contourint(1000,abcissa(i),fun);
    Wrap1time(i) = toc();
    tic();
    Wrap2Estimate(i) = contourint2(-0.1,0+abcissa(i),-0.1,0.1,.0001,fun);
    Wrap2time(i) = toc();
    fun = 'dml16(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissa(i),0.5);
    Weekstime(i) = toc();
end

error = True-Wrap1Estimate;
Wrap1Error(8,:) = [sum(Wrap1time),max(Wrap1time),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-Wrap2Estimate;
Wrap2Error(8,:) = [sum(Wrap2time),max(Wrap2time),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-BromEstimate;
BromError(8,:) = [sum(Bromtime),max(Bromtime),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-WeeksEstimate;
WeeksError(8,:) = [sum(Weekstime),max(Weekstime),mean(abs(error)),norm(error,2),max(abs(error))];

figure
plot(T,True,'LineWidth',4,'Color','black')
hold on
%plot(T,real(Wrap1Estimate),'o','MarkerSize',7,'Color','green')
plot(T,real(Wrap2Estimate),'o','MarkerSize',7,'Color','cyan')
plot(T,real(BromEstimate),'o','MarkerSize',7,'Color','red')
plot(T,real(WeeksEstimate),'o','MarkerSize',7,'Color','blue')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','square trapez','Bromwich adapt','Weeks')
hold off

% A&V fcn 5

True = av5(T);

for i=1:100
    fun = @(x)avl5(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*10,abcissa(i)+1i*10);
    Bromtime(i) = toc();
    tic();
    Wrap1Estimate(i) = contourint(1000,abcissa(i),fun);
    Wrap1time(i) = toc();
    tic();
    Wrap2Estimate(i) = contourint2(-0.1,0+abcissa(i),-0.1,0.1,.0001,fun);
    Wrap2time(i) = toc();
    fun = 'avl5(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissa(i),0.5);
    Weekstime(i) = toc();
end

error = True-Wrap1Estimate;
Wrap1Error(9,:) = [sum(Wrap1time),max(Wrap1time),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-Wrap2Estimate;
Wrap2Error(9,:) = [sum(Wrap2time),max(Wrap2time),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-BromEstimate;
BromError(9,:) = [sum(Bromtime),max(Bromtime),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-WeeksEstimate;
WeeksError(9,:) = [sum(Weekstime),max(Weekstime),mean(abs(error)),norm(error,2),max(abs(error))];

figure
plot(T,True,'LineWidth',4,'Color','black')
hold on
plot(T,real(Wrap1Estimate),'o','MarkerSize',7,'Color','green')
plot(T,real(Wrap2Estimate),'o','MarkerSize',7,'Color','cyan')
plot(T,real(BromEstimate),'o','MarkerSize',7,'Color','red')
plot(T,real(WeeksEstimate),'o','MarkerSize',7,'Color','blue')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','circular trapez','square trapez','Bromwich adapt','Weeks')
hold off

% A&V fcn 10

True = av10(T);

for i=1:100
    fun = @(x)avl10(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*10,abcissa(i)+1i*10);
    Bromtime(i) = toc();
    tic();
    Wrap1Estimate(i) = contourint(1000,abcissa(i),fun);
    Wrap1time(i) = toc();
    tic();
    Wrap2Estimate(i) = contourint2(-0.1,0+abcissa(i),-0.1,0.1,.0001,fun);
    Wrap2time(i) = toc();
    fun = 'avl10(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissa(i),0.5);
    Weekstime(i) = toc();
end

error = True-Wrap1Estimate;
Wrap1Error(10,:) = [sum(Wrap1time),max(Wrap1time),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-Wrap2Estimate;
Wrap2Error(10,:) = [sum(Wrap2time),max(Wrap2time),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-BromEstimate;
BromError(10,:) = [sum(Bromtime),max(Bromtime),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-WeeksEstimate;
WeeksError(10,:) = [sum(Weekstime),max(Weekstime),mean(abs(error)),norm(error,2),max(abs(error))];

figure
plot(T,True,'LineWidth',4,'Color','black')
hold on
plot(T,real(Wrap1Estimate),'o','MarkerSize',7,'Color','green')
plot(T,real(Wrap2Estimate),'o','MarkerSize',7,'Color','cyan')
plot(T,real(BromEstimate),'o','MarkerSize',7,'Color','red')
plot(T,real(WeeksEstimate),'o','MarkerSize',7,'Color','blue')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','circular trapez','square trapez','Bromwich adapt','Weeks')
hold off

% A&V fcn 12

True = av12(T);

for i=1:100
    fun = @(x)avl12(x).*exp(T(i).*x)./(2*pi*1i);
    tic();
    BromEstimate(i) = integral(fun,abcissa(i)-1i*10,abcissa(i)+1i*10);
    Bromtime(i) = toc();
    tic();
    Wrap1Estimate(i) = contourint(1000,abcissa(i),fun);
    Wrap1time(i) = toc();
    tic();
    Wrap2Estimate(i) = contourint2(-0.1,0+abcissa(i),-0.1,0.1,.0001,fun);
    Wrap2time(i) = toc();
    fun = 'avl12(s)';
    tic();
    WeeksEstimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissa(i),0.5);
    Weekstime(i) = toc();
end

error = True-Wrap1Estimate;
Wrap1Error(11,:) = [sum(Wrap1time),max(Wrap1time),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-Wrap2Estimate;
Wrap2Error(11,:) = [sum(Wrap2time),max(Wrap2time),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-BromEstimate;
BromError(11,:) = [sum(Bromtime),max(Bromtime),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-WeeksEstimate;
WeeksError(11,:) = [sum(Weekstime),max(Weekstime),mean(abs(error)),norm(error,2),max(abs(error))];

figure
plot(T,True,'LineWidth',4,'Color','black')
hold on
plot(T,real(Wrap1Estimate),'o','MarkerSize',7,'Color','green')
plot(T,real(Wrap2Estimate),'o','MarkerSize',7,'Color','cyan')
plot(T,real(BromEstimate),'o','MarkerSize',7,'Color','red')
plot(T,real(WeeksEstimate),'o','MarkerSize',7,'Color','blue')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','circular trapez','square trapez','Bromwich adapt','Weeks')
hold off

format long g

Wrap1Error
Wrap2Error
BromError
WeeksError
