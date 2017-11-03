function F = avl10(S)
%s=0
[n,m]=size(S);
if(n>m)
    sz=n;
else
    sz=m;
end
F = zeros(n,m);
for i=1:sz
F(i) = exp(S(i))*besselk(1,S(i))/S(i);
end

