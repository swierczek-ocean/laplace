function f = ilapt10(T,a)
[n,m]=size(T);
if(n>m)
    sz=n;
else
    sz=m;
end
f = zeros(n,m);
for i=1:sz
f(i) = (1/6)*T(i)*T(i)*T(i)*exp(a*T(i));
end

