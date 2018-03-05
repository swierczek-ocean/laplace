function f = ilapt72(T,a)
[n,m] = size(T);
f = zeros(n,m);
Tfloor = floor(T);
for ii=1:n
    for jj=1:m
        for kk=1:Tfloor(ii,jj)
            f(ii,jj) = f(ii,jj) + a^kk;
        end
    end
end
end


