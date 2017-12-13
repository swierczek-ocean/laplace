function f = ilapt29(T)
[n,m]=size(T);
if(n>m)
    sz=n;
else
    sz=m;
end
f = zeros(n,m);
for i=1:sz
f(i) = (3*T(i)*sin(T(i))+(T(i)^2)*cos(T(i)))/8;
end

