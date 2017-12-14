
tic()
set(0,'DefaultFigureVisible','on')
colors
ep = 2;
cof = 0.5;
bromerr = [];
weekserr = [];
b=0.8;
T=2:2:200;
num = max(size(T));
sigma = cof./T.^ep;
sigmaW = sigma;
B = 75;
sw=0;
BromEstimate = zeros(1,num);
WeeksEstimate = zeros(1,num);


% D&M fcn 1

True = dm1(T);

for i=1:num
    fun = 'dml1(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)dml1(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('DM 1')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
print('DM1','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('DM 1 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
print('DM1error','-djpeg')
hold off


% D&M fcn 2

True = dm2(T);

for i=1:num
    fun = 'dml2(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)dml2(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('DM 2')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('DM2','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('DM 2 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('DM2error','-djpeg')
hold off

% D&M fcn 5

True = dm5(T);

for i=1:num
    fun = 'dml5(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)dml5(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('DM 5')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('DM5','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('DM 5 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('DM5error','-djpeg')
hold off

% D&M fcn 6

True = dm6(T);

for i=1:num
    fun = 'dml6(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)dml6(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('DM 6')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('DM6','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('DM 6 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('DM6error','-djpeg')
hold off


% D&M fcn 8

True = dm8(T);

for i=1:num
    fun = 'dml8(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)dml8(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('DM 8')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('DM8','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('DM 8 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('DM8error','-djpeg')
hold off

% D&M fcn 9

True = dm9(T);

for i=1:num
    fun = 'dml9(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)dml9(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('DM 9')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('DM9','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('DM 9 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('DM9error','-djpeg')
hold off

% D&M fcn 10

True = dm10(T);

for i=1:num
    fun = 'dml10(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)dml10(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('DM 10')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('DM10','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('DM 10 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('DM10error','-djpeg')
hold off

% D&M fcn 11

True = dm11(T);

for i=1:num
    fun = 'dml11(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)dml11(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('DM 11')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('DM11','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('DM 11 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('DM11error','-djpeg')
hold off

close all

% D&M fcn 13

True = dm13(T);

for i=1:num
    fun = 'dml13(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)dml13(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('DM 13')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('DM13','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('DM 13 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('DM13error','-djpeg')
hold off

% D&M fcn 15

True = dm15(T);

for i=1:num
    fun = 'dml15(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)dml15(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('DM 15')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('DM15','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('DM 15 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('DM15error','-djpeg')
hold off

% D&M fcn 16

True = dm16(T);

for i=1:num
    fun = 'dml16(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)dml16(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('DM 16')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('DM16','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('DM 16 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('DM16error','-djpeg')
hold off

% A&V fcn 5

True = av5(T);

for i=1:num
    fun = 'avl5(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)avl5(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('AV 5')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('AV5','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('AV 5 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('AV5error','-djpeg')
hold off

% A&V fcn 10

True = av10(T);

for i=1:num
    fun = 'avl10(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)avl10(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('AV 10')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('AV10','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('AV 10 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('AV10error','-djpeg')
hold off

close all

% S fcn 1

True = ilapt1(T);

for i=1:num
    fun = 'lapt1(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)lapt1(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 1')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('ILT1','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 1 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('ILT1error','-djpeg')
hold off


% S fcn 2

True = ilapt2(T);

for i=1:num
    fun = 'lapt2(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)lapt2(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 2')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('ILT2','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 2 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('ILT2error','-djpeg')
hold off


% S fcn 3

True = ilapt3(T);

for i=1:num
    fun = 'lapt3(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)lapt3(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 3')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('ILT3','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 3 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('ILT3error','-djpeg')
hold off


% S fcn 4

True = ilapt4(T);

for i=1:num
    fun = 'lapt4(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)lapt4(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 4')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('ILT4','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 4 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('ILT4error','-djpeg')
hold off


% S fcn 5

True = ilapt5(T);

for i=1:num
    fun = 'lapt5(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)lapt5(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 5')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('ILT5','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 5 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('ILT5error','-djpeg')
hold off


% S fcn 12

True = ilapt12(T);

for i=1:num
    fun = 'lapt12(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)lapt12(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 12')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('ILT12','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 12 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('ILT12error','-djpeg')
hold off

close all

% S fcn 13

True = ilapt13(T);

for i=1:num
    fun = 'lapt13(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)lapt13(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 13')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('ILT13','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 13 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('ILT13error','-djpeg')
hold off


% S fcn 14

True = ilapt14(T);

for i=1:num
    fun = 'lapt14(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)lapt14(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 14')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('ILT14','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 14 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('ILT14error','-djpeg')
hold off


% S fcn 15

True = ilapt15(T);

for i=1:num
    fun = 'lapt15(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)lapt15(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 15')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('ILT15','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 15 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('ILT5error','-djpeg')
hold off


% S fcn 16

True = ilapt16(T);

for i=1:num
    fun = 'lapt16(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)lapt16(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 16')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('ILT16','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 16 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('ILT16error','-djpeg')
hold off

% S fcn 17

True = ilapt17(T);

for i=1:num
    fun = 'lapt17(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)lapt17(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 17')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('ILT17','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 17 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('ILT17error','-djpeg')
hold off

close all

% S fcn 18

True = ilapt18(T);

for i=1:num
    fun = 'lapt18(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)lapt18(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 18')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('ILT18','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 18 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('ILT18error','-djpeg')
hold off


% S fcn 19

True = ilapt19(T);

for i=1:num
    fun = 'lapt19(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)lapt19(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 19')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('ILT19','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 19 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('ILT19error','-djpeg')
hold off


% S fcn 20

True = ilapt20(T);

for i=1:num
    fun = 'lapt20(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)lapt20(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 20')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('ILT20','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 20 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('ILT20error','-djpeg')
hold off


% S fcn 26

True = ilapt26(T);

for i=1:num
    fun = 'lapt26(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)lapt26(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 26')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('ILT26','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 26 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('ILT26error','-djpeg')
hold off

% S fcn 27

True = ilapt27(T);

for i=1:num
    fun = 'lapt27(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)lapt27(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 27')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('ILT27','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 27 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('ILT27error','-djpeg')
hold off

% S fcn 28

True = ilapt28(T);

for i=1:num
    fun = 'lapt28(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)lapt28(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 28')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('ILT28','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 28 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('ILT28error','-djpeg')
hold off


% S fcn 29

True = ilapt29(T);

for i=1:num
    fun = 'lapt29(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)lapt29(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 29')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('ILT29','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 29 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('ILT29error','-djpeg')
hold off


% S fcn 30

True = ilapt30(T);

for i=1:num
    fun = 'lapt30(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)lapt30(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 30')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('ILT30','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 30 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('ILT30error','-djpeg')
hold off


% S fcn 31

True = ilapt31(T);

for i=1:num
    fun = 'lapt31(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)lapt31(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 31')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('ILT31','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 31 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('ILT31error','-djpeg')
hold off

% S fcn 32

True = ilapt32(T);

for i=1:num
    fun = 'lapt32(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)lapt32(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 32')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('ILT32','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 32 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('ILT32error','-djpeg')
hold off

% S fcn 33

True = ilapt33(T);

for i=1:num
    fun = 'lapt33(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)lapt33(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 33')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('ILT33','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 33 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('ILT33error','-djpeg')
hold off


% S fcn 34

True = ilapt34(T);

for i=1:num
    fun = 'lapt34(s)';
    [WeeksEstimate(i),OutParamS,OutErrorS,LaguCoeff] = WeeksMethod(fun,T(i),0.0000001,sw);
    OutParamS = OutParamS
    
    
    fun = @(y)lapt34(sigma(i)+1i*y).*fmint(y,sigma(i),T(i));
    BromEstimate(i) = integral(fun,-B,B,'AbsTol',1e-16);
end

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,real(BromEstimate),'+','MarkerSize',8,'Color',Color(:,9))
plot(T,real(WeeksEstimate),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 34')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt','Weeks')
%fullscreen
print('ILT34','-djpeg')
hold off

WeeksRelError = abs(WeeksEstimate-True)./abs(True);
BromRelError = abs(BromEstimate-True)./abs(True);
bromerr = [bromerr,mean(BromRelError)];
weekserr = [weekserr,mean(WeeksRelError)];

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
hold on
semilogy(T,real(WeeksRelError),'*','MarkerSize',8,'Color',Color(:,12))
title('ILT 34 Error')
xlabel('time')
ylabel('f(t)')
legend('Bromwich error','Weeks error')
%fullscreen
print('ILT34error','-djpeg')
hold off



close all


bromerr = bromerr
weekserr = weekserr

toc()