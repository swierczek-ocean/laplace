function f = dm8(T)
[n,m]=size(T);
if(n>m)
    sz=n;
else
    sz=m;
end
f = zeros(n,m);
for i=1:sz
f(i) = sin(T(i));
end

