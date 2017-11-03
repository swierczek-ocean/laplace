function f = dml2(T)
[n,m]=size(T);
if(n>m)
    sz=n;
else
    sz=m;
end
f = zeros(n,m);
for i=1:sz
f(i) = (1/sqrt(pi*T(i)))*cos(2*sqrt(T(i)));
end
