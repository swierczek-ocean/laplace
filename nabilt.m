function [f] = nabilt(fun,t,ub,jj,sw,Sings)

%% Naive Adaptive Bromwich Contour
% calculates inverse Laplace transform by numerical integration
% along a finite portion of the Browmich contour
% after translating the function to have its rightmost
% singularity to have real part zero

%% find the singularities
if(sw==1)
    [hshift,vshift] = singfind(fun);
else
    hshift = Sings(jj,1);
    vshift = Sings(jj,2);
    %[hshift,vshift] = Sings(jj,:);
end
%%

%% test that prior part worked
fprintf('horizontal shift = %g for function %g\n',hshift,jj)
fprintf('vertical shift = %g for function %g\n',vshift,jj)
%%


%% elementary transformation
fun2 = @(x)fun(x+hshift);

%% loop calculation
[n,m] = size(t);
f = zeros(n,m);

for ii=1:n
    for jj=1:m
        T = t(ii,jj);
        sigma = 0.5/T^2;
        fun3 = @(x)(fun2(sigma + 1i.*x).*iltint(x,sigma,T));
        fun4 = @(x)double(fun3(x));
        f(ii,jj) = real(exp(sigma*T)*integral(fun4,-ub+vshift,ub+vshift,'AbsTol',1e-13)/2/pi);
    end
end

f = double(exp(hshift.*t).*f);

end

