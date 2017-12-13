function f = ilapt8(T,a)
[n,m]=size(T);
if(n>m)
    sz=n;
else
    sz=m;
end
f = zeros(n,m);
for i=1:sz
f(i) = T(i)*exp(a*T(i));
end

