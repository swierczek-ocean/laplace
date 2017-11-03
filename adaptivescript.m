T=1:1:20;

True = dm6(T);
abcissa = 10;

Estimate = zeros(1,20);


for i=1:20
    fun = @(x)dml6(x).*exp(T(i).*x);
    %Estimate(i) = integral(fun,abcissa-1i*50,abcissa+1i*50);
    %Estimate(i) = contourint(100,abcissa,fun);
    Estimate(i) = contourint2(-1,1,-1,1,.0005,fun);
end

error = True-Estimate;
figure
plot(T,True,'*','MarkerSize',5,'Color','blue')
hold on
plot(T,real(Estimate),'o','MarkerSize',5,'Color','red')


figure
plot(T,abs(error),'o','MarkerSize',5,'Color','red')

Estimate'