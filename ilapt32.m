function f = ilapt32(T)
[n,m]=size(T);
if(n>m)
    sz=n;
else
    sz=m;
end
f = zeros(n,m);
for i=1:sz
f(i) = 0.5*(T(i)^2)*sin(T(i));
end


