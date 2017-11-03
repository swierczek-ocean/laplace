function F = dml4(S)
%s=-0.2
[n,m]=size(S);
if(n>m)
    sz=n;
else
    sz=m;
end
F = zeros(n,m);
for i=1:sz
F(i) = 1/((S(i)+0.2)*(S(i)+0.2)+1);
end

