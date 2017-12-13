function F = lapt11(s)
% sing @ s = a
[n,m]=size(s);
if(n>m)
    sz=n;
else
    sz=m;
end
F = zeros(n,m);
for i=1:sz
F(i) = 1/(s(i)-1)^5;
end

