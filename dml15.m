function F = dml15(S)
%s=0?
[n,m]=size(S);
if(n>m)
    sz=n;
else
    sz=m;
end
F = zeros(n,m);
for i=1:sz
F(i) = exp(-4*sqrt(S(i)));
end

