function F = dml3(S)
%s=-0.5
[n,m]=size(S);
if(n>m)
    sz=n;
else
    sz=m;
end
F = zeros(n,m);
for i=1:sz
F(i) = 1/(S(i)+0.5);
end

