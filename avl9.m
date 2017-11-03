function F = avl9(S)
%s=0
[n,m]=size(S);
if(n>m)
    sz=n;
else
    sz=m;
end
F = zeros(n,m);
for i=1:sz
F(i) = log(S(i)-1+sqrt(S(i)*S(i)-2*S(i)))/sqrt(S(i)*S(i)-2*S(i));
end

