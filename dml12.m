function F = dml12(S)
%s=0
[n,m]=size(S);
if(n>m)
    sz=n;
else
    sz=m;
end
F = zeros(n,m);
for i=1:sz
F(i) = (1/S(i))*tanh(S(i)/2);
end

