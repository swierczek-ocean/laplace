function [f] = nabilt(fun,t)

%% Naive Adaptive Bromwich Contour
% calculates inverse Laplace transform by numerical integration
% along a finite portion of the Browmich contour
% after translating the function to have its rightmost
% singularity to have real part zero
tic()

%% find the singularities
syms s
P = poles(fun,s);
fun1 = @(s)1/fun(s);
Z = vpasolve(fun1(s) == 0, s);
S = [P;Z];
S = real(S);
S = sort(S,'descend');
shift = S(1)
%%

%% elementary transformation
fun2 = @(s)fun(s+shift);


%% loop calculation
[n,m] = size(t);
f = zeros(n,m);
ub = 80;

for ii=1:n
    for jj=1:m
        T = t(ii,jj);
        sigma = 0.5/T^2;
        fun3 = @(x)fun2(sigma + 1i*x)*iltint(x,sigma,T);
        f(ii,jj) = exp(sigma*T)*integral(fun3,-ub,ub)/2/pi/1i;
    end
end

f = exp(shift.*t).*f;

toc()
end

