function F = lapt35(s)
% sing @ s = +-i1
[n,m]=size(s);
if(n>m)
    sz=n;
else
    sz=m;
end
F = zeros(n,m);
for i=1:sz
F(i) = (s(i)^4-6*(s(i)^2)+1)/(s(i)^2+1)^4;
end

