function [f] = nabilt(fun,t,ub)

%% Naive Adaptive Bromwich Contour
% calculates inverse Laplace transform by numerical integration
% along a finite portion of the Browmich contour
% after translating the function to have its rightmost
% singularity to have real part zero

%% find the singularities
syms s
P = poles(fun,s);
fun1 = @(s)1/fun(s);
Z = vpasolve(fun1(s) == 0, s);
S = [P;Z];
Q = real(S);
[Q,I] = sort(Q,'descend');
if(isempty(Q)==1)
    hshift = 0;
    vshift = 0;
else
    hshift = Q(1);
    vshift = double(imag(S(I(1))));
end
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

