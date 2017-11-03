function f = dm14(T)
[n,m]=size(T);
if(n>m)
    sz=n;
else
    sz=m;
end
f = zeros(n,m);
for i=1:sz
f(i) = (exp(-0.25*T(i))-exp(-0.5*T(i)))/(sqrt(4*pi*T(i)*T(i)*T(i)));
end

