function f = ilapt200(T,a,n)
[r,m]=size(T);
if(r>m)
    sz=r;
else
    sz=m;
end
f = zeros(r,m);
for i=1:sz
f(i) = (T(i)^(n-1))*exp(a*T(i))/gamma(n);
end

