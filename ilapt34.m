function f = ilapt34(T,a)
f = ((3+a^2.*T.^2).*sinh(a.*T)+5*a.*T.*cosh(a.*T))./(8*a);
end

