T=1:1:20;
abcissa = 1./T;
Estimate = zeros(1,20);



True = dm6(T);

for i=1:20
    x = T(i);
    fun = 'dml6(s)';
    %Estimate(i) = integral(fun,abcissa(i)-1i*50,abcissa(i)+1i*50);
    %Estimate(i) = contourint(500,abcissa(i),fun);
    %Estimate(i) = contourint2(-1,1,-1,1,.0005,fun);
    Estimate(i) = wfnWeeksCoreSigmab(fun,T(i),256,abcissa(i),0.5);
end

error = True-Estimate;
figure
plot(T,True,'*','MarkerSize',5,'Color','blue')
hold on
plot(T,real(Estimate),'o','MarkerSize',5,'Color','red')


figure
plot(T,abs(error),'o','MarkerSize',5,'Color','red')

Estimate'