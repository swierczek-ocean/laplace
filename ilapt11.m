function f = ilapt11(T)
[n,m]=size(T);
if(n>m)
    sz=n;
else
    sz=m;
end
f = zeros(n,m);
for i=1:sz
f(i) = (1/24)*(T(i)^4)*exp(T(i));
end

