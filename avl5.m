function F = avl5(S)
%s=0
[n,m]=size(S);
if(n>m)
    sz=n;
else
    sz=m;
end
F = zeros(n,m);
for i=1:sz
F(i) = (exp(-0.25/S(i)))/(S(i)^(3/2));
end

