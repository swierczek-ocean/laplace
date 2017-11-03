function F = avl1(S)
%s=?
[n,m]=size(S);
if(n>m)
    sz=n;
else
    sz=m;
end
F = zeros(n,m);
for i=1:sz
F(i) = 1/(sqrt(S(i))+sqrt(S(i)+1));
end

