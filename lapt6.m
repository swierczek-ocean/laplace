function F = lapt6(s)
% sing @ s = 1
[n,m]=size(s);
if(n>m)
    sz=n;
else
    sz=m;
end
F = zeros(n,m);
for i=1:sz
F(i) = 1/(s(i)-1);
end

