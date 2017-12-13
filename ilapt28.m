function f = ilapt28(T)
[n,m]=size(T);
if(n>m)
    sz=n;
else
    sz=m;
end
f = zeros(n,m);
for i=1:sz
f(i) = ((1+(T(i)^2))*sin(T(i))-T(i)*cos(T(i)))/8;
end

