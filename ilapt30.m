function f = ilapt30(T)
[n,m]=size(T);
if(n>m)
    sz=n;
else
    sz=m;
end
f = zeros(n,m);
for i=1:sz
f(i) = ((3-(T(i)^2))*sin(T(i))+5*T(i)*cos(T(i)))/8;
end

