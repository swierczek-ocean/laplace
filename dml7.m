function F = dml7(S)
%s=-1
[n,m]=size(S);
if(n>m)
    sz=n;
else
    sz=m;
end
F = zeros(n,m);
for i=1:sz
F(i) = 1/(S(i)+1)^2;
end

