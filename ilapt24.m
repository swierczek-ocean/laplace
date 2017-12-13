function f = ilapt24(T)
[n,m]=size(T);
if(n>m)
    sz=n;
else
    sz=m;
end
f = zeros(n,m);
for i=1:sz
f(i) = cosh(T(i))+0.5*T(i)*sinh(T(i));
end

