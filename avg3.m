function f = avg3(T)
[n,m]=size(T);
if(n>m)
    sz=n;
else
    sz=m;
end
f = zeros(n,m);
for i=1:sz
f(i) = sin(T(i))/sqrt(2*pi*T(i)*T(i)*T(i));
end

