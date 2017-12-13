function F = lapt4(s)
% sing @ s = 0
[n,m]=size(s);
if(n>m)
    sz=n;
else
    sz=m;
end
F = zeros(n,m);
for i=1:sz
F(i) = 1/s(i)^4;
end

