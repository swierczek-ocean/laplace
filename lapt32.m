function F = lapt32(s)
% sing @ s = +-ia
[n,m]=size(s);
if(n>m)
    sz=n;
else
    sz=m;
end
F = zeros(n,m);
for i=1:sz
F(i) = (3*(s(i)^2)-1)/(s(i)^2+1)^3;
end
