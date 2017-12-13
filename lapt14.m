function F = lapt14(s)
% sing @ s = a,-a
[n,m]=size(s);
if(n>m)
    sz=n;
else
    sz=m;
end
F = zeros(n,m);
for i=1:sz
F(i) = 1/(s(i)^2-1);
end

