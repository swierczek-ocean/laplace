T=1:1:100;
abcissa = 1./T;

Wrap1Error = zeros(14,5);
Wrap2Error = zeros(14,5);
BromError = zeros(14,5);
WeeksError = zeros(14,5);

Wrap1Estimate = zeros(1,100);
Wrap2Estimate = zeros(1,100);
BromEstimate = zeros(1,100);
WeeksEstimate = zeros(1,100);

Wrap1time = zeros(1,100);
Wrap2time = zeros(1,100);
Bromtime = zeros(1,100);
Weekstime = zeros(1,100);

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
end

error = True-Wrap1Estimate;
Wrap1Error(1,:) = [sum(Wrap1time),max(Wrap1time),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-Wrap2Estimate;
Wrap2Error(1,:) = [sum(Wrap2time),max(Wrap2time),mean(abs(error)),norm(error,2),max(abs(error))];

error = True-BromEstimate;
BromError(1,:) = [sum(Bromtime),max(Bromtime),mean(abs(error)),norm(error,2),max(abs(error))];







