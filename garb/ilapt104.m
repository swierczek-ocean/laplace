function f = ilapt4(T)
[n,m]=size(T);
if(n>m)
    sz=n;
else
    sz=m;
end
f = zeros(n,m);
for i=1:sz
f(i) = (1/6)*T(i)*T(i)*T(i);
end

