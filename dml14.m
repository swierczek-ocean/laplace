function F = dml14(S)
[n,m]=size(S);
if(n>m)
    sz=n;
else
    sz=m;
end
F = zeros(n,m);
for i=1:sz
F(i) = sqrt(S(i)+0.5) - sqrt(S(i)+0.25);
end

