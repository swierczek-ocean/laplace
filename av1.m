function f = av1(T)
[n,m]=size(T);
if(n>m)
    sz=n;
else
    sz=m;
end
f = zeros(n,m);
for i=1:sz
f(i) = (1-exp(-T(i)))/(2*sqrt(pi*T(i)*T(i)*T(i)));
end

