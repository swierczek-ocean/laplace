function f = av5(T)
[n,m]=size(T);
if(n>m)
    sz=n;
else
    sz=m;
end
f = zeros(n,m);
for i=1:sz
f(i) = 2*sin(sqrt(T(i)))/sqrt(pi);
end

