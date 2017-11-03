function f = avg5(T)
[n,m]=size(T);
if(n>m)
    sz=n;
else
    sz=m;
end
f = zeros(n,m);
for i=1:sz
f(i) = besselj(1,sqrt(T(i)*T(i)+2*T(i)))/sqrt(T(i)*T(i)+2*T(i));
end

