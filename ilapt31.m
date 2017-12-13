function f = ilapt31(T)
[n,m]=size(T);
if(n>m)
    sz=n;
else
    sz=m;
end
f = zeros(n,m);
for i=1:sz
f(i) = ((8-(T(i)^2))*cos(T(i))-7*T(i)*sin(T(i)))/8;
end

