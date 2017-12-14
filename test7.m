
tic()
set(0,'DefaultFigureVisible','on')
colors
ep = 2;
cof = 0.5;
brommean = [];
weeksmean = [];
b=0.8;
T=51.3:10:301.3;
num = max(size(T));
sigma = cof./T.^ep;
sigmaW = sigma;
B = 75;
BromError = [];
WeeksError = [];
BromEstimate = zeros(1,num);
WeeksEstimate = zeros(1,num);
Bromtime = zeros(1,num);
Weekstime = zeros(1,num);


% D&M fcn 1

True = dm1(T);

for i=1:num
    fun = 'dml1(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,0);
    OutParamS = OutParamS
    
    
    fun = @(y)dml1(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
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
plot(T,True,'LineWidth',3,'Color',Color(:,28))
hold on
plot(T,real(BromEstimate),'o','MarkerSize',7,'Color',Color(:,22))
plot(T,real(WeeksEstimate),'*','MarkerSize',7,'Color',Color(:,28))
title('DM 1')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);

figure
semilogy(T,real(BromRelError),'o','MarkerSize',7,'Color',Color(:,22))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',7,'Color',Color(:,28))
title('DM 1 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
hold off

