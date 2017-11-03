function F = avl3(S)
%s=1
[n,m]=size(S);
if(n>m)
    sz=n;
else
    sz=m;
end
F = zeros(n,m);
for i=1:sz
F(i) = (S(i)-sqrt(S(i)*S(i)-1))/sqrt(S(i)*S(i)-1);
end

