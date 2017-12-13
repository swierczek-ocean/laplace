function f = ilapt17(T)
[n,m]=size(T);
if(n>m)
    sz=n;
else
    sz=m;
end
f = zeros(n,m);
for i=1:sz
f(i) = 0.5*T(i)*sin(T(i));
end

