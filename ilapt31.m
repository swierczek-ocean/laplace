function f = ilapt31(T,a)
f = ((a.*T.^2).*cosh(a.*T)-T.*sinh(a.*T))./(8*a^3);
end

